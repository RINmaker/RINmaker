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

#include "bonds.h"
#include "chemical_entity.h"

#include "rin_params.h"
#include "log_manager.h"
#include "spatial/kdtree.h"

namespace fs = std::filesystem;

using lm = log_manager;

using chemical_entity::aminoacid, chemical_entity::component,
        chemical_entity::atom, chemical_entity::ring, chemical_entity::ionic_group;

using std::vector, std::string, std::list, std::set, std::map, std::unordered_map, std::function, std::optional,
        std::shared_ptr, std::make_shared, std::make_unique, std::is_base_of, std::ifstream, std::runtime_error,
        std::pair, std::nullopt;

using rin::parameters, prelude::interval;

struct rin::maker::impl
{
public:
    vector<aminoacid const*> _aminoacids;

    kdtree<atom, 3> _hdonor_tree, _vdw_tree;
    vector<atom const*> _hacceptor_vector, _vdw_vector, _cation_vector;

    kdtree<ring, 3> _ring_tree, _pication_ring_tree;
    vector<ring const*> _ring_vector, _pication_ring_vector;

    kdtree<ionic_group, 3> _positive_ion_tree;
    vector<ionic_group const*> _negative_ion_vector;

    kdtree<atom, 3> _alpha_carbon_tree, _beta_carbon_tree;
    vector<atom const*> _alpha_carbon_vector, _beta_carbon_vector;

    // ss bonds are directly parsed, not computed by us
    vector<shared_ptr<bond::ss const>> _ss_bonds;

    string _pdb_name;

    ~impl()
    { for (auto* res: _aminoacids) delete res; }
};

vector<shared_ptr<rin::maker>> rin::maker::parse_models(fs::path const& pdb_path)
{
    ifstream pdb_file;
    pdb_file.open(pdb_path);

    // might throw
    if (!pdb_file.is_open())
        throw runtime_error("could not open " + pdb_path.string() + "\n");

    string ls;
    uint32_t ln = 0;
    vector<pair<uint32_t, string>> numbered_lines;
    while (getline(pdb_file, ls))
        numbered_lines.emplace_back(ln++, ls);
    pdb_file.close();

    // atom records are split in MODEL sections
    vector<vector<records::atom>> models;

    // other records are in global sections
    vector<records::ss> ssbond_records;
    vector<records::helix> helix_records;
    vector<records::sheet_piece> sheet_records;

    vector<records::atom> tmp_model;
    for (auto const& line: numbered_lines)
    {
        auto const& line_str = line.second;
        auto const& line_num = line.first;

        auto const record_type = prelude::trim(line_str.substr(0, 6));

        // TODO exception
        // MODEL should always balance ENDMDL
        // and there cannot be 2 consecutive MODEL or 2 consecutive ENDMDL
        // but there can be 0 MODEL and 0 ENDMDL

        if (record_type == "ATOM")
            tmp_model.emplace_back(line_str, line_num);

        else if (record_type == "HELIX")
            helix_records.emplace_back(line_str, line_num);

        else if (record_type == "SHEET")
            sheet_records.emplace_back(line_str, line_num);

        else if (record_type == "SSBOND")
            ssbond_records.emplace_back(line_str, line_num);

        else if (record_type == "ENDMDL")
            models.push_back(tmp_model);

        // else if (record_type == "MODEL");
    }

    // TODO this can only happen when no MODEL/ENDMDL is specified
    if (!tmp_model.empty())
        models.push_back(tmp_model);

    auto const pdb_name = pdb_path.stem().string();
    vector<shared_ptr<rin::maker>> rin_makers;
    for (auto const& model: models)
        rin_makers.emplace_back(
                make_shared<rin::maker>(pdb_name, model, ssbond_records, helix_records, sheet_records));

    return rin_makers;
}

template<typename Record>
class secondary_structure_helper final
{
private:
    unordered_map<string, map<interval<int>, Record, interval<int>::less>> _map;

public:
    void insert(Record const& record)
    {
        auto chain = _map.find(record.init_chain_id());
        if (chain == _map.end())
        {
            map<interval<int>, Record, interval<int>::less> new_chain;
            new_chain.insert({record.range(), record});
            _map.insert({record.init_chain_id(), new_chain});
        }
        else
        { chain->second.insert({record.range(), record}); }
    }

    [[nodiscard]]
    size_t size() const
    { return _map.size(); }

    optional<Record> find(aminoacid const& res) const
    {
        auto chain = _map.find(res.chain_id());
        if (chain != _map.end())
        {
            auto kv = chain->second.find(interval<int>(res.sequence_number(), res.sequence_number()));
            if (kv != chain->second.end())
                return kv->second;
        }

        return nullopt;
    }
};

rin::maker::maker(string const& pdb_name,
                  vector<records::atom> const& atom_records,
                  vector<records::ss> const& ssbond_records,
                  vector<records::helix> const& helix_records,
                  vector<records::sheet_piece> const& sheet_records) : pimpl{make_unique<impl>()}
{
    vector<records::atom> tmp_atoms;
    vector<aminoacid*> tmp_aminoacids;

    lm::main()->info("parsing pdb lines...");

    for (auto const& record: atom_records)
    {
        if (!tmp_atoms.empty() && !record.same_res(tmp_atoms.back()))
        {
            tmp_aminoacids.emplace_back(new aminoacid(tmp_atoms, pdb_name));
            tmp_atoms.clear();
        }

        tmp_atoms.push_back(record);
    }

    if (!tmp_atoms.empty())
        tmp_aminoacids.emplace_back(new aminoacid(tmp_atoms, pdb_name));

    for (auto const& record: ssbond_records)
        pimpl->_ss_bonds.emplace_back(make_shared<bond::ss>(record));

    auto sheet_helper = secondary_structure_helper<records::sheet_piece>();
    for (auto const& record: sheet_records)
        sheet_helper.insert(record);

    auto helix_helper = secondary_structure_helper<records::helix>();
    for (auto const& record: helix_records)
        helix_helper.insert(record);

    lm::main()->info("finding the appropriate secondary structure for each aminoacid...");

    if (sheet_helper.size() != 0 || helix_helper.size() != 0)
    {
        for (auto res: tmp_aminoacids)
        {
            res->make_secondary_structure();

            auto sheet = sheet_helper.find(*res);
            if (sheet.has_value())
                res->make_secondary_structure(sheet.value());

            auto helix = helix_helper.find(*res);
            if (helix.has_value())
                res->make_secondary_structure(helix.value());
        }
    }

    pimpl->_aminoacids.reserve(tmp_aminoacids.size());
    for (auto res: tmp_aminoacids)
        pimpl->_aminoacids.emplace_back(res);

    lm::main()->info("retrieving components from aminoacids...");

    // used only to build the corresponding kdtrees
    vector<atom const*> hdonors;
    vector<ionic_group const*> positives;

    for (auto const* res: pimpl->_aminoacids)
    {

        auto carbon = res->ca();
        if (carbon != nullptr)
            pimpl->_alpha_carbon_vector.push_back(carbon);

        carbon = res->cb();
        if (carbon != nullptr)
            pimpl->_beta_carbon_vector.push_back(carbon);

        for (auto const* a: res->atoms())
        {
            if (a->is_a_hydrogen_donor())
                hdonors.push_back(a);

            if (a->is_a_hydrogen_acceptor())
                pimpl->_hacceptor_vector.push_back(a);

            if (a->is_a_vdw_candidate())
                pimpl->_vdw_vector.push_back(a);

            if (a->is_a_cation())
                pimpl->_cation_vector.push_back(a);
        }

        auto const* group = res->positive_ionic_group();
        if (group != nullptr)
            positives.push_back(group);

        group = res->negative_ionic_group();
        if (group != nullptr)
            pimpl->_negative_ion_vector.push_back(group);

        auto const* ring = res->primary_ring();
        if (ring != nullptr)
        {
            pimpl->_ring_vector.push_back(ring);

            if (ring->is_a_pication_candidate())
                pimpl->_pication_ring_vector.push_back(ring);
        }

        ring = res->secondary_ring();
        if (ring != nullptr)
        {
            pimpl->_ring_vector.push_back(ring);

            if (ring->is_a_pication_candidate())
                pimpl->_pication_ring_vector.push_back(ring);
        }
    }

    lm::main()->info("hydrogen acceptors: {}", pimpl->_hacceptor_vector.size());
    lm::main()->info("hydrogen donors: {}", hdonors.size());
    lm::main()->info("vdw candidates: {}", pimpl->_vdw_vector.size());
    lm::main()->info("cations: {}", pimpl->_cation_vector.size());
    lm::main()->info("aromatic rings: {}", pimpl->_ring_vector.size());
    lm::main()->info("aromatic rings (valid for pication): {}", pimpl->_ring_vector.size());

    lm::main()->info("building acceleration structures...");

    pimpl->_hdonor_tree = kdtree<atom, 3>(hdonors);
    pimpl->_vdw_tree = kdtree<atom, 3>(pimpl->_vdw_vector);

    pimpl->_ring_tree = kdtree<ring, 3>(pimpl->_ring_vector);
    pimpl->_pication_ring_tree = kdtree<ring, 3>(pimpl->_pication_ring_vector);
    pimpl->_positive_ion_tree = kdtree<ionic_group, 3>(positives);

    pimpl->_alpha_carbon_tree = kdtree<atom, 3>(pimpl->_alpha_carbon_vector);
    pimpl->_beta_carbon_tree = kdtree<atom, 3>(pimpl->_beta_carbon_vector);
}

template<typename Bond, typename Entity1, typename Entity2>
vector<shared_ptr<Bond const>>
find_bonds(vector<Entity1 const*> const& vec, kdtree<Entity2, 3> const& tree, double dist, parameters const& params)
{
    static_assert(
            is_base_of<component, Entity1>::value,
            "template typename Entity1 must inherit from type chemical_entity::component");
    static_assert(
            is_base_of<component, Entity2>::value,
            "template typename Entity2 must inherit from type chemical_entity::component");
    static_assert(
            is_base_of<bond::base, Bond>::value, "template typename BondFunc must inherit from type bond::base");

    vector<shared_ptr<Bond const>> bonds;
    for (auto e1: vec)
    {
        auto neighbors = tree.range_search(*e1, dist);
        for (auto e2: neighbors)
        {
            auto bond = Bond::test(params, *e1, *e2);
            if (bond != nullptr)
                bonds.emplace_back(bond);
        }
    }

    return bonds;
}

vector<shared_ptr<bond::base const>> filter_hbond_realistic(vector<shared_ptr<bond::base const>> const& input)
{
    std::set<shared_ptr<bond::hydrogen const>> hydrogen_bonds_output;
    std::vector<shared_ptr<bond::hydrogen const>> hydrogen_bonds_input;
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
    auto can_be_added = [&](shared_ptr<bond::hydrogen const> bond) -> bool
    {
        return (get_bond_count(donors_bond_count, bond->donor_ptr()) <
                bond->donor().how_many_hydrogen_can_donate() &&
                get_bond_count(hydrogen_bond_count, bond->hydrogen_ptr()) < 1 &&
                //An hydrogen can make only one bond
                get_bond_count(acceptors_bond_count, bond->acceptor_ptr()) <
                bond->acceptor().how_many_hydrogen_can_accept());
    };
    auto add_bond = [&](shared_ptr<bond::hydrogen const> bond) -> void
    {
        inc_bond_count(donors_bond_count, bond->donor_ptr());
        inc_bond_count(hydrogen_bond_count, bond->hydrogen_ptr());
        inc_bond_count(acceptors_bond_count, bond->acceptor_ptr());
        hydrogen_bonds_output.insert(bond);
    };

    //Extract hydrogen bonds from input
    for (auto& i: input)
    {
        if (i->get_type() == "hydrogen")
            hydrogen_bonds_input.push_back(std::dynamic_pointer_cast<bond::hydrogen const>(i));
    }

    //Order from smallest to largest energy
    sort(
            hydrogen_bonds_input.begin(), hydrogen_bonds_input.end(),
            [](shared_ptr<bond::hydrogen const> a, shared_ptr<bond::hydrogen const> b)
            { return a->get_energy() < b->get_energy(); });

    //Add as many hydrogen bonds as possible
    for (auto i: hydrogen_bonds_input)
    {
        if (can_be_added(i))
            add_bond(i);
    }

    //Let's build the output list
    vector<shared_ptr<bond::base const>> output;
    for (auto i: input)
    {
        //Insert i into the output if it is not an hydrogen or if it is in the filtered list
        if (i->get_type() != "hydrogen" ||
            hydrogen_bonds_output.find(std::dynamic_pointer_cast<bond::hydrogen const>(i)) !=
            hydrogen_bonds_output.end())
        {
            output.push_back(i);
        }
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
        auto const hydrogen_bonds = find_bonds<bond::hydrogen>(
                pimpl->_hacceptor_vector,
                pimpl->_hdonor_tree,
                params.query_dist_hbond(),
                params);

        auto const vdw_bonds = remove_duplicates(
                find_bonds<bond::vdw>(
                        pimpl->_vdw_vector,
                        pimpl->_vdw_tree,
                        params.query_dist_vdw(),
                        params));

        auto const ionic_bonds = find_bonds<bond::ionic>(
                pimpl->_negative_ion_vector,
                pimpl->_positive_ion_tree,
                params.query_dist_ionic(),
                params);

        auto const cationpi_bonds = find_bonds<bond::pication>(
                pimpl->_cation_vector,
                pimpl->_pication_ring_tree,
                params.query_dist_pica(),
                params);

        auto const pipistack_bonds = remove_duplicates(
                find_bonds<bond::pipistack>(
                        pimpl->_ring_vector,
                        pimpl->_ring_tree,
                        params.query_dist_pipi(),
                        params));

        switch (params.network_policy())
        {
        case parameters::network_policy_t::ALL:
            append(results, hydrogen_bonds);
            append(results, vdw_bonds);
            append(results, ionic_bonds);
            append(results, cationpi_bonds);
            append(results, pipistack_bonds);
            append(results, pimpl->_ss_bonds);
            break;

        case parameters::network_policy_t::BEST_ONE:
            append(results, hydrogen_bonds);
            append(results, vdw_bonds);
            append(results, ionic_bonds);
            append(results, cationpi_bonds);
            append(results, pipistack_bonds);
            append(results, pimpl->_ss_bonds);

            results = filter_best(results);
            break;

        case parameters::network_policy_t::BEST_PER_TYPE:
            append(results, filter_best(hydrogen_bonds));
            append(results, filter_best(vdw_bonds));
            append(results, filter_best(ionic_bonds));
            append(results, filter_best(cationpi_bonds));
            append(results, filter_best(pipistack_bonds));
            append(results, filter_best(pimpl->_ss_bonds));
            break;
        }

        if (params.hbond_realistic())
            results = filter_hbond_realistic(results);

        break;
    }

    case parameters::interaction_type_t::GENERIC_ALPHA:
    {
        auto alpha_bonds = find_bonds<bond::generic_bond>(
                pimpl->_alpha_carbon_vector,
                pimpl->_alpha_carbon_tree,
                params.query_dist_alpha(),
                params);

        append(results, filter_best(alpha_bonds));
        break;
    }

    case parameters::interaction_type_t::GENERIC_BETA:
    {
        auto beta_bonds = find_bonds<bond::generic_bond>(
                pimpl->_beta_carbon_vector,
                pimpl->_beta_carbon_tree,
                params.query_dist_beta(),
                params);

        append(results, filter_best(beta_bonds));
        break;
    }
    }

    // results = remove_duplicates(results);
    lm::main()->info("there are {} valid bonds after filtering", results.size());

    return {pimpl->_pdb_name, params, pimpl->_aminoacids, results};
}
