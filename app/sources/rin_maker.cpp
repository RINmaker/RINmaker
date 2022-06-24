#include "rin_maker.h"

#include <string>
#include <list>
#include <map>
#include <unordered_map>

#include <functional>
#include <optional>

#include <fstream>

#include <exception>
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

using rin::parameters, prelude::interval;


/*rin::maker::maker(string const& pdb_name,
                  vector<record::atom> const& atom_records,
                  vector<record::ss> const& ssbond_records,
                  vector<record::helix> const& helix_records,
                  vector<record::sheet_piece> const& sheet_records)*/
rin::maker::maker(gemmi::Model const& model, gemmi::Structure const& structure)
{
    // we are filling the private implementation piece-by-piece, so we need a non-const temporary here
    // at the end of the constructor we will store it in the private member pimpl, which is a const*
    auto tmp_pimpl = make_shared<rin::maker::impl>();

    // GEMMI already parses all records and then groups atoms together in the respective residues
    // so all information is already here. We will just need to rebuild rings and ionic groups
    for (auto const& chain : model.chains)
        for (auto const& residue : chain.residues)
            tmp_pimpl->aminoacids.emplace_back(residue, chain, model, structure);

    lm::main()->info("retrieving components from aminoacids...");

    // these are used only to build the corresponding kdtrees
    vector<atom> hdonors;
    vector<ionic_group> positives;

    for (auto const& res: tmp_pimpl->aminoacids)
    {

        auto const ca = res.get_alpha_carbon();
        if (ca.has_value())
            tmp_pimpl->alpha_carbon_vector.push_back(*ca);

        auto const cb = res.get_beta_carbon();
        if (cb.has_value())
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

        auto const pos_group = res.get_positive_ionic_group();
        if (pos_group.has_value())
            positives.push_back(*pos_group);

        auto const neg_group = res.get_negative_ionic_group();
        if (neg_group.has_value())
            tmp_pimpl->negative_ion_vector.push_back(*neg_group);

        auto const ring_setup = [&](std::optional<ring> const& ring)
        {
            if (ring.has_value())
            {
                tmp_pimpl->ring_vector.push_back(*ring);

                if (ring->is_a_pication_candidate())
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
    lm::main()->info("aromatic rings: {}", tmp_pimpl->ring_vector.size());
    lm::main()->info("aromatic rings (valid for pication): {}", tmp_pimpl->ring_vector.size());

    lm::main()->info("building acceleration structures...");

    tmp_pimpl->hdonor_tree = kdtree<atom, 3>(hdonors);
    tmp_pimpl->vdw_tree = kdtree<atom, 3>(tmp_pimpl->vdw_vector);

    tmp_pimpl->ring_tree = kdtree<ring, 3>(tmp_pimpl->ring_vector);
    tmp_pimpl->pication_ring_tree = kdtree<ring, 3>(tmp_pimpl->pication_ring_vector);
    tmp_pimpl->positive_ion_tree = kdtree<ionic_group, 3>(positives);

    tmp_pimpl->alpha_carbon_tree = kdtree<atom, 3>(tmp_pimpl->alpha_carbon_vector);
    tmp_pimpl->beta_carbon_tree = kdtree<atom, 3>(tmp_pimpl->beta_carbon_vector);

    for (auto const& connection : structure.connections)
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
            is_base_of<aminoacid::component, Entity1>::value,
            "template typename Entity1 must inherit from type chemical_entity::aminoacid::component");
    static_assert(
            is_base_of<aminoacid::component, Entity2>::value,
            "template typename Entity2 must inherit from type chemical_entity::aminoacid::component");
    static_assert(
            is_base_of<bond::base, Bond>::value, "template typename Bond must inherit from type bond::base");

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
        return (get_bond_count(donors_bond_count, &bond->donor()) <
                bond->donor().how_many_hydrogen_can_donate() &&
                get_bond_count(hydrogen_bond_count, &bond->hydrogen_atom()) < 1 &&
                //An hydrogen can make only one bond
                get_bond_count(acceptors_bond_count, &bond->acceptor()) <
                bond->acceptor().how_many_hydrogen_can_accept());
    };
    auto add_bond = [&](shared_ptr<bond::hydrogen const> const& bond) -> void
    {
        inc_bond_count(donors_bond_count, &bond->donor());
        inc_bond_count(hydrogen_bond_count, &bond->hydrogen_atom());
        inc_bond_count(acceptors_bond_count, &bond->acceptor());
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
        auto const pair_id = b->get_id_simple();
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

        break;
    }

    case parameters::interaction_type_t::CONTACT_MAP:
    {
        vector<shared_ptr<bond::generic_bond const>> generic_bonds;
        switch (params.cmap_type())
        {
        case rin::parameters::contact_map_type_t::ALPHA:
            generic_bonds = find_bonds<bond::generic_bond>(
                    pimpl->alpha_carbon_vector,
                    pimpl->alpha_carbon_tree,
                    params.query_dist_cmap(),
                    params);
            break;

        case rin::parameters::contact_map_type_t::BETA:
            generic_bonds = find_bonds<bond::generic_bond>(
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

    // results = remove_duplicates(results);
    lm::main()->info("there are {} valid bonds after filtering", results.size());

    return {pimpl->pdb_name, params, pimpl->aminoacids, results};
}
