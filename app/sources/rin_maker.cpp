#include "rin_maker.h"

#include <string>
#include <list>
#include <map>
#include <set>
#include <unordered_map>

#include <optional>

#include <utility>

#include "ns_bond.h"
#include "ns_chemical_entity.h"

#include "rin_params.h"
#include "log_manager.h"
#include "spatial/kdtree.h"

#include "private/impl_rin_maker.h"

namespace fs = std::filesystem;

using lm = log_manager;

using chemical_entity::aminoacid, chemical_entity::atom, chemical_entity::ring, chemical_entity::ionic_group;

using std::vector, std::string, std::list, std::set, std::map, std::unordered_map, std::function, std::optional,
        std::shared_ptr, std::make_shared, std::make_unique, std::is_base_of, std::ifstream, std::runtime_error,
        std::pair, std::nullopt;

using rin::parameters;

/**
 * To be used in sorted structures, such as std::map (see std::map requirements @ cppreference.com).
 * <br/>
 * <br/>
 * Semantics:
 * <br/>
 * - a < b means "interval a comes before interval b and they do not overlap"
 * <br/>
 * - a > b means b < a
 * <br/>
 * - a == b means overlapping.
 */
template<typename T>
struct interval final
{
private:
    T _inf;
    T _sup;

public:
    interval(T const &a, T const &b) : _inf(a < b ? a : b), _sup(a < b ? b : a)
    {}

    struct less final
    {
        bool operator()(interval const &lhs, interval const &rhs) const
        { return lhs._sup < rhs._inf; }
    };
};

template <typename Secondary>
struct secondary_structure_helper_map
{
private:
    unordered_map<string, map<interval<int>, Secondary, interval<int>::less>> chain_to_map;

    static std::optional<std::tuple<std::string, interval<int>>> maybe_chain_and_interval(Secondary sstructure, gemmi::Model const& model)
    {
        auto cra_start = model.find_cra(sstructure.start);
        auto cra_end = model.find_cra(sstructure.end);

        if (cra_start.residue != nullptr && cra_end.residue != nullptr)
        {
            auto seq_start = cra_start.residue->seqid.num.value;
            auto seq_end = cra_end.residue->seqid.num.value;

            return {{cra_start.chain != nullptr ? cra_start.chain->name : "", {seq_start, seq_end}}};
        }
        else
        { return std::nullopt; }
    }

public:
    void maybe_insert(Secondary const& t, gemmi::Model const& model)
    {
        auto maybe_stuff = maybe_chain_and_interval(t, model);

        if (!maybe_stuff.has_value())
            return;

        auto chain_name = std::get<0>(*maybe_stuff);
        auto sequence_interval = std::get<1>(*maybe_stuff);

        auto chain_and_ts = chain_to_map.find(chain_name);
        if (chain_and_ts == chain_to_map.end())
        {
            map<interval<int>, Secondary, interval<int>::less> tmp{{sequence_interval, t}};
            chain_to_map.insert({chain_name, tmp});
        }
        else
        { chain_and_ts->second.insert_or_assign(sequence_interval, t); }
    }

    std::optional<Secondary> maybe_find(gemmi::Residue const& residue, gemmi::Chain const& chain)
    {
        std::optional<Secondary> result{};
        if (auto chain_and_ts = chain_to_map.find(chain.name); chain_and_ts != chain_to_map.end())
        {
            if (auto t = chain_and_ts->second.find(interval(residue.seqid.num.value, residue.seqid.num.value));
                t != chain_and_ts->second.end())
                result = {t->second};
        }
        return result;
    }

    auto empty() const
    { return chain_to_map.empty(); }
};

/**
 * This function tries to check that all parsed residues as actual aminoacids.
 * <br/>
 * For now it just issues a warning if it does not find any.
 * @param residues
 */
void warn_if_not_protein(std::vector<chemical_entity::aminoacid> const& residues)
{
    static std::set<std::string, std::less<>> names{
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        "SEC", "PYL"
    };

    for(auto const& res : residues)
        if(auto const match = names.find(res.get_name()); match != names.end())
            return;

    lm::main()->warn("it seems that the pdb/cif file does not contain any valid aminoacid!");
}

rin::maker::maker(gemmi::Model const& model, gemmi::Structure const& protein,  rin::parameters const& params)
{
    secondary_structure_helper_map<gemmi::Helix> helix_map;
    for (auto const& helix: protein.helices)
        helix_map.maybe_insert(helix, model);

    secondary_structure_helper_map<gemmi::Sheet::Strand> strand_map;
    for (auto const& sheet: protein.sheets)
        for (auto const& strand: sheet.strands)
            strand_map.maybe_insert(strand, model);

    // we are filling the private implementation piece-by-piece, so we need a non-const temporary here
    // at the end of the constructor we will store it in the private member pimpl, which is a const*
    auto tmp_pimpl = make_shared<rin::maker::impl>();

    auto try_build_aminoacid =
        [&helix_map, &strand_map, &tmp_pimpl, &model, &protein, &params]
        (auto const& residue, auto const& chain)
    {
        try
        {
            if (helix_map.empty() && strand_map.empty())
                tmp_pimpl->aminoacids.emplace_back(residue, chain, model, protein, params);
            else
            {
                std::optional<std::variant<gemmi::Helix, gemmi::Sheet::Strand>> maybe_sstruct{std::nullopt};
                if (!(maybe_sstruct = helix_map.maybe_find(residue, chain)).has_value())
                    maybe_sstruct = strand_map.maybe_find(residue, chain);

                tmp_pimpl->aminoacids.emplace_back(residue, chain, model, protein, params, maybe_sstruct);
            }
        }
        catch (std::exception const& e)
        {
            if (params.illformed_policy() == rin::parameters::illformed_policy_t::FAIL)
            {
                lm::main()->error("aborting: {}", e.what());
                throw;
            }
            else
                lm::main()->warn("skipping residue: {}", e.what());
        }
    };

    lm::main()->info("building aminoacids...");

    // GEMMI already parses all records and then groups atoms together in the respective residues
    // so all information is already here. We will just need to rebuild rings and ionic groups
    for (auto const& chain: model.chains)
        for (auto const& residue: chain.residues)
            if (!residue.is_water() || !params.skip_water())
                try_build_aminoacid(residue, chain);

    warn_if_not_protein(tmp_pimpl->aminoacids);

    lm::main()->info("extracting ionic groups, rings and other entities...");

    // these are used only to build the corresponding kdtrees
    vector<atom> hdonors;
    vector<ionic_group> positives;

    for (auto const& res: tmp_pimpl->aminoacids)
    {
        if (auto const& ca = res.get_alpha_carbon(); ca.has_value())
            tmp_pimpl->alpha_carbon_vector.push_back(*ca);

        if (auto const& cb = res.get_beta_carbon(); cb.has_value())
            tmp_pimpl->beta_carbon_vector.push_back(*cb);

        for (auto const& a : res.get_atoms())
        {
            if (a.is_hydrogen_donor())
                hdonors.push_back(a);

            if (a.is_hydrogen_acceptor())
                tmp_pimpl->hacceptor_vector.push_back(a);

            if (a.is_vdw_candidate())
                tmp_pimpl->vdw_vector.push_back(a);

            if (a.is_cation())
                tmp_pimpl->cation_vector.push_back(a);
        }

        if (auto const& pos_group = res.get_positive_ionic_group(); pos_group.has_value())
            positives.push_back(*pos_group);

        if (auto const& neg_group = res.get_negative_ionic_group(); neg_group.has_value())
            tmp_pimpl->negative_ion_vector.push_back(*neg_group);

        auto const ring_setup = [&](std::optional<ring> const& ring)
        {
            if (ring.has_value())
            {
                tmp_pimpl->ring_vector.push_back(*ring);

                if (ring->is_pication_candidate())
                    tmp_pimpl->pication_ring_vector.push_back(*ring);
            }
        };

        ring_setup(res.get_primary_ring());
        ring_setup(res.get_secondary_ring());
    }

    lm::main()->info("hydrogen acceptors: {}", tmp_pimpl->hacceptor_vector.size());
    lm::main()->info("hydrogen donors: {}", hdonors.size());
    lm::main()->info("vdw candidates: {}", tmp_pimpl->vdw_vector.size());
    lm::main()->info("cations: {}", tmp_pimpl->cation_vector.size());
    lm::main()->info("aromatic rings (total): {}", tmp_pimpl->ring_vector.size());
    lm::main()->info("aromatic rings (cation-pi only): {}", tmp_pimpl->pication_ring_vector.size());

    lm::main()->info("building kdtrees...");

    tmp_pimpl->hdonor_tree = kdtree<atom, 3>(hdonors);
    tmp_pimpl->vdw_tree = kdtree<atom, 3>(tmp_pimpl->vdw_vector);

    tmp_pimpl->ring_tree = kdtree<ring, 3>(tmp_pimpl->ring_vector);
    tmp_pimpl->pication_ring_tree = kdtree<ring, 3>(tmp_pimpl->pication_ring_vector);
    tmp_pimpl->positive_ion_tree = kdtree<ionic_group, 3>(positives);

    tmp_pimpl->alpha_carbon_tree = kdtree<atom, 3>(tmp_pimpl->alpha_carbon_vector);
    tmp_pimpl->beta_carbon_tree = kdtree<atom, 3>(tmp_pimpl->beta_carbon_vector);

    for (auto const& connection : protein.connections)
        if (connection.type == gemmi::Connection::Type::Disulf)
            tmp_pimpl->ss_bonds.push_back(std::make_shared<bond::ss>(connection));

    pimpl = tmp_pimpl;
}

rin::maker::~maker() = default;

template<typename Bond, typename Entity1, typename Entity2>
vector<shared_ptr<Bond const>>
find_bonds(vector<Entity1> const& vec, kdtree<Entity2, 3> const& tree, double dist, parameters const& params)
{
    static_assert(
        std::is_base_of_v<aminoacid::component, Entity1>,
        "template typename Entity1 must inherit from type chemical_entity::aminoacid::component");
    static_assert(
        std::is_base_of_v<aminoacid::component, Entity2>,
        "template typename Entity2 must inherit from type chemical_entity::aminoacid::component");
    static_assert(
        std::is_base_of_v<bond::base, Bond>,
        "template typename Bond must inherit from type bond::base");

    vector<shared_ptr<Bond const>> bonds;
    for (auto const& e1 : vec)
    {
        auto neighbors = tree.range_search(e1, dist);
        for (auto const& e2 : neighbors)
        {
            auto bond = Bond::test(params, e1, e2);
            if (bond != nullptr)
                bonds.emplace_back(bond);
        }
    }

    return bonds;
}

std::vector<shared_ptr<bond::hydrogen const>> filter_hbond_realistic(std::vector<shared_ptr<bond::hydrogen const>> input)
{
    std::vector<shared_ptr<bond::hydrogen const>> output;
    std::unordered_map<chemical_entity::atom const*, int> donors_bond_count;
    std::unordered_map<chemical_entity::atom const*, int> hydrogen_bond_count;
    std::unordered_map<chemical_entity::atom const*, int> acceptors_bond_count;

    //Get the bonds count of an atom
    auto get_bond_count = [](
            std::unordered_map<chemical_entity::atom const*, int>& container, chemical_entity::atom const* atom) -> int
    {
        if (container.find(atom) == container.end())
            return 0;
        else
            return container[atom];
    };

    //Increase the bonds count of an atom
    auto inc_bond_count = [](
            std::unordered_map<chemical_entity::atom const*, int>& container,
            chemical_entity::atom const* atom) -> void
    {
        if (container.find(atom) == container.end())
            container[atom] = 0;
        container[atom]++;
    };
    auto can_be_added = [&](shared_ptr<bond::hydrogen const> const& bond) -> bool
    {
        return (get_bond_count(donors_bond_count, &bond->get_donor()) <
                bond->get_donor().how_many_hydrogen_can_donate() &&
                get_bond_count(hydrogen_bond_count, &bond->get_hydrogen_atom()) < 1 &&
                //An hydrogen can make only one bond
                get_bond_count(acceptors_bond_count, &bond->get_acceptor()) <
                bond->get_acceptor().how_many_hydrogen_can_accept());
    };
    auto add_bond = [&](shared_ptr<bond::hydrogen const> const& bond) -> void
    {
        inc_bond_count(donors_bond_count, &bond->get_donor());
        inc_bond_count(hydrogen_bond_count, &bond->get_hydrogen_atom());
        inc_bond_count(acceptors_bond_count, &bond->get_acceptor());
        output.push_back(bond);
    };

    //Order from smallest to largest energy
    sort(input.begin(), input.end(),
         [](shared_ptr<bond::hydrogen const> const& a, shared_ptr<bond::hydrogen const> const& b)
         { return a->get_energy() < b->get_energy(); });

    //Add as many hydrogen bonds as possible
    for (const auto& i: input)
    {
        if (can_be_added(i))
            add_bond(i);
    }

    return output;
}

template<typename Bond>
vector<shared_ptr<Bond const>> remove_duplicates(vector<shared_ptr<Bond const>> const& unfiltered)
{
    vector<shared_ptr<Bond const>> results;
    results.reserve(unfiltered.size());

    set<string> unique_ids;
    for (auto const& b: unfiltered)
    {
        auto const bond_id = b->get_id();
        if (unique_ids.find(bond_id) == unique_ids.end())
        {
            unique_ids.insert(bond_id);
            results.push_back(b);
        }
    }

    return results;
}

template<typename Bond>
vector<shared_ptr<Bond const>> filter_best(vector<shared_ptr<Bond const>> const& unfiltered)
{
    vector<shared_ptr<Bond const>> results;
    results.reserve(unfiltered.size());

    unordered_map<string, shared_ptr<Bond const>> res_pairs;
    for (auto const& b: unfiltered)
    {
        auto const pair_id = b->get_pair_id();
        auto const pair = res_pairs.find(pair_id);
        if (pair == res_pairs.end() || *b < *pair->second)
            res_pairs.insert_or_assign(pair_id, b);
    }

    for (auto const& pair: res_pairs)
        results.push_back(pair.second);

    return results;
}

template<typename Bond>
void append(vector<shared_ptr<bond::base const>>& dst, vector<shared_ptr<Bond const>> const& src)
{ dst.insert(dst.end(), src.begin(), src.end()); }

rin::graph rin::maker::operator()(parameters const& params) const
{
    vector<shared_ptr<bond::base const>> results;
    switch (params.interaction_type())
    {
    case parameters::interaction_type_t::NONCOVALENT_BONDS:
    {
        lm::main()->info("finding all bonds...");
        auto hydrogen_bonds = find_bonds<bond::hydrogen>(
                pimpl->hacceptor_vector,
                pimpl->hdonor_tree,
                params.query_dist_hbond(),
                params);
        if (params.hbond_realistic())
            hydrogen_bonds = filter_hbond_realistic(hydrogen_bonds);

        auto const vdw_bonds = remove_duplicates(
                find_bonds<bond::vdw>(
                        pimpl->vdw_vector,
                        pimpl->vdw_tree,
                        params.query_dist_vdw(),
                        params));

        auto const ionic_bonds = find_bonds<bond::ionic>(
                pimpl->negative_ion_vector,
                pimpl->positive_ion_tree,
                params.query_dist_ionic(),
                params);

        auto const pication_bonds = find_bonds<bond::pication>(
                pimpl->cation_vector,
                pimpl->pication_ring_tree,
                params.query_dist_pica(),
                params);

        auto const pipistack_bonds = remove_duplicates(
                find_bonds<bond::pipistack>(
                        pimpl->ring_vector,
                        pimpl->ring_tree,
                        params.query_dist_pipi(),
                        params));

        switch (params.network_policy())
        {
        case parameters::network_policy_t::ALL:
            append(results, hydrogen_bonds);
            append(results, vdw_bonds);
            append(results, ionic_bonds);
            append(results, pication_bonds);
            append(results, pipistack_bonds);
            append(results, pimpl->ss_bonds);
            break;

        case parameters::network_policy_t::BEST_ONE:
            append(results, hydrogen_bonds);
            append(results, vdw_bonds);
            append(results, ionic_bonds);
            append(results, pication_bonds);
            append(results, pipistack_bonds);
            append(results, pimpl->ss_bonds);

            results = filter_best(results);
            break;

        case parameters::network_policy_t::BEST_PER_TYPE:
            append(results, filter_best(hydrogen_bonds));
            append(results, filter_best(vdw_bonds));
            append(results, filter_best(ionic_bonds));
            append(results, filter_best(pication_bonds));
            append(results, filter_best(pipistack_bonds));
            append(results, filter_best(pimpl->ss_bonds));
            break;
        }

        // hydrophobic bonds are just put into the rin _after_ fltering
        append(results, remove_duplicates(find_bonds<bond::hydrophobic>(
            pimpl->alpha_carbon_vector,
            pimpl->alpha_carbon_tree,
            7.5,
            params
        )));
        break;
    }

    case parameters::interaction_type_t::CONTACT_MAP:
    {
        lm::main()->info("generating contact map...");
        vector<shared_ptr<bond::contact const>> generic_bonds{};
        switch (params.cmap_type())
        {
        case rin::parameters::contact_map_type_t::ALPHA:
            generic_bonds = find_bonds<bond::contact>(
                    pimpl->alpha_carbon_vector,
                    pimpl->alpha_carbon_tree,
                    params.query_dist_cmap(),
                    params);
            break;

        case rin::parameters::contact_map_type_t::BETA:
            generic_bonds = find_bonds<bond::contact>(
                    pimpl->beta_carbon_vector,
                    pimpl->beta_carbon_tree,
                    params.query_dist_cmap(),
                    params);
            break;
        }

        append(results, filter_best(generic_bonds));
        break;
    }
    }

    lm::main()->info("count: {}", results.size());

    return {pimpl->pdb_name, pimpl->aminoacids, results};
}
