#include "rin_maker.h"

#include <functional>
#include <fstream>
#include <exception>

#include <string>
#include <list>

#include <optional>

#include <map>
#include <unordered_map>

#include "rin_network.h"
#include "rin_params.h"
#include "bonds.h"

using std::list;
using std::string;
using std::unordered_map;
using std::map;
using std::function;

template<typename Record>
class secondary_structure_helper final
{
private:
    std::unordered_map<string, std::map<interval<int>, Record, interval<int>::less>> _map;

public:
    void insert(Record const& record)
    {
        auto chain = _map.find(record.init_chain_id());
        if (chain == _map.end())
        {
            std::map<interval<int>, Record, interval<int>::less> new_chain;
            new_chain.insert({record.range(), record});
            _map.insert({record.init_chain_id(), new_chain});
        }
        else
        {
            chain->second.insert({record.range(), record});
        }
    }

    [[nodiscard]]
    int size() const
    { return _map.size(); }

    std::optional<Record> find(chemical_entity::aminoacid& res) const
    {
        auto chain = _map.find(res.chain_id());
        if (chain != _map.end())
        {
            auto kv = chain->second.find(interval<int>(res.sequence_number(), res.sequence_number()));
            if (kv != chain->second.end())
                return kv->second;
        }

        return std::nullopt;
    }
};

rin::maker::maker(std::string const& pdb_name, std::vector<numbered_line_t>::iterator const begin, std::vector<numbered_line_t>::iterator const end) : _pdb_name(pdb_name)
{
    auto sheet_records = secondary_structure_helper<records::sheet_piece>();
    auto helix_records = secondary_structure_helper<records::helix>();

    vector<records::atom> tmp_atoms;

    lm::main()->info("parsing pdb lines...");

    for (auto it = begin; it != end; ++it)
    {
        auto line = it->second;
        auto line_number = it->first;

        auto const record_type = prelude::trim(line.substr(0, 6));
        if (record_type == "ATOM")
        {
            // atoms are grouped by aminoacid: we can simply fill a collection until we change residue
            records::atom record(line, line_number);
            if (!tmp_atoms.empty() && !record.same_res(tmp_atoms.back()))
            {
                _aminoacids.push_back(new aminoacid(tmp_atoms, pdb_name));
                tmp_atoms.clear();
            }

            tmp_atoms.push_back(record);
        }

        else if (record_type == "HELIX")
            helix_records.insert(records::helix(line, line_number));

        else if (record_type == "SHEET")
            sheet_records.insert(records::sheet_piece(line, line_number));

        else if (record_type == "SSBOND")
            _ss_bonds.emplace_back(new bond::ss(records::ss(line, line_number)));
    }

    if (!tmp_atoms.empty())
        _aminoacids.push_back(new aminoacid(tmp_atoms, pdb_name));

    lm::main()->info("finding the appropriate secondary structure for each aminoacid...");

    if (sheet_records.size() != 0 || helix_records.size() != 0)
    {
        for (auto& res: _aminoacids)
        {
            res->make_secondary_structure();

            auto sheet = sheet_records.find(*res);
            if (sheet.has_value())
                res->make_secondary_structure(sheet.value());

            auto helix = helix_records.find(*res);
            if(helix.has_value())
                res->make_secondary_structure(helix.value());
        }
    }

    lm::main()->info("retrieving components from aminoacids...");

    // used only to build the corresponding kdtrees
    vector<atom const*> hdonors;
    vector<ionic_group const*> positives;

    for (auto const* res: _aminoacids)
    {

        auto carbon = res->ca();
        if (carbon != nullptr)
            _alpha_carbon_vector.push_back(carbon);

        carbon = res->cb();
        if (carbon != nullptr)
            _beta_carbon_vector.push_back(carbon);

        for (auto const* a: res->atoms())
        {
            if (a->is_a_hydrogen_donor())
                hdonors.push_back(a);

            if (a->is_a_hydrogen_acceptor())
                _hacceptor_vector.push_back(a);

            if (a->is_a_vdw_candidate())
                _vdw_vector.push_back(a);

            if (a->is_a_cation())
                _cation_vector.push_back(a);
        }

        auto const* group = res->positive_ionic_group();
        if (group != nullptr)
            positives.push_back(group);

        group = res->negative_ionic_group();
        if (group != nullptr)
            _negative_ion_vector.push_back(group);

        auto const* ring = res->primary_ring();
        if (ring != nullptr)
        {
            _ring_vector.push_back(ring);

            if (ring->is_a_pication_candidate())
                _pication_ring_vector.push_back(ring);
        }

        ring = res->secondary_ring();
        if (ring != nullptr)
        {
            _ring_vector.push_back(ring);

            if (ring->is_a_pication_candidate())
                _pication_ring_vector.push_back(ring);
        }
    }

    lm::main()->info("hydrogen acceptors: {}", _hacceptor_vector.size());
    lm::main()->info("hydrogen donors: {}", hdonors.size());
    lm::main()->info("vdw candidates: {}", _vdw_vector.size());
    lm::main()->info("cations: {}", _cation_vector.size());
    lm::main()->info("aromatic rings: {}", _ring_vector.size());
    lm::main()->info("aromatic rings (valid for pication): {}", _ring_vector.size());

    lm::main()->info("building acceleration structures...");

    _hdonor_tree = kdtree<chemical_entity::atom, 3>(hdonors);
    _vdw_tree = kdtree<chemical_entity::atom, 3>(_vdw_vector);

    _ring_tree = kdtree<chemical_entity::ring, 3>(_ring_vector);
    _pication_ring_tree = kdtree<chemical_entity::ring, 3>(_pication_ring_vector);
    _positive_ion_tree = kdtree<chemical_entity::ionic_group, 3>(positives);

    _alpha_carbon_tree = kdtree<atom, 3>(_alpha_carbon_vector);
    _beta_carbon_tree = kdtree<atom, 3>(_beta_carbon_vector);
}

rin::maker::~maker()
{ for (auto* res: _aminoacids) delete res; }

using chemical_entity::component;

template<typename Bond, typename Entity1, typename Entity2>
static void find_bonds(
        network& net, std::vector<Entity1 const*> const& vec, kdtree<Entity2, 3> const& tree, double dist, rin::parameters const& params)
{
    static_assert(
            std::is_base_of<component, Entity1>::value, "template typename Entity1 must inherit from type entity::component");
    static_assert(
            std::is_base_of<component, Entity2>::value, "template typename Entity2 must inherit from type entity::component");
    static_assert(
            std::is_base_of<bond::computed, Bond>::value, "template typename BondFunc must inherit from type bond::computed");

    for (auto* e1: vec)
    {
        auto neighbors = tree.range_search(*e1, dist);
        for (auto* e2: neighbors)
            Bond::test(net, params, *e1, *e2);
    }
}

rin::graph rin::maker::operator()(parameters const& params) const
{
    network _network;

    list<bond::base const*> results;
    switch (params.interaction_type())
    {
    case parameters::interaction_type_t::NONCOVALENT_BONDS:
        find_bonds<bond::hydrogen>(
                _network, _hacceptor_vector, _hdonor_tree, params.query_dist_hbond(), params);

        find_bonds<bond::vdw>(
                _network, _vdw_vector, _vdw_tree, params.query_dist_vdw(), params);

        find_bonds<bond::ionic>(
                _network, _negative_ion_vector, _positive_ion_tree, params.query_dist_ionic(), params);

        find_bonds<bond::pication>(
                _network, _cation_vector, _pication_ring_tree, params.query_dist_pica(), params);

        find_bonds<bond::pipistack>(
                _network, _ring_vector, _ring_tree, params.query_dist_pipi(), params);

        switch (params.network_policy())
        {
        case parameters::network_policy_t::ALL:
            results = _network.get_all();
            break;

        case parameters::network_policy_t::BEST_PER_TYPE:
            results = _network.get_multiple();
            break;

        case parameters::network_policy_t::BEST_ONE:
            results = _network.get_one();
            break;
        }

        // TODO
        // think of a way to place the filter before get_multiple/get_all...
        // (maybe it is necessary to create a new net copying only filtered results).
        if (params.hbond_realistic())
            results = _network.filter_hbond_realistic(results);
        break;

    case parameters::interaction_type_t::ALPHA_BACKBONE:
        find_bonds<bond::generico>(
                _network, _alpha_carbon_vector, _alpha_carbon_tree, params.query_dist_alpha(), params);
        results = _network.get_one();
        break;

    case parameters::interaction_type_t::BETA_BACKBONE:
        find_bonds<bond::generico>(
                _network, _beta_carbon_vector, _beta_carbon_tree, params.query_dist_beta(), params);
        results = _network.get_one();
        break;
    }

    lm::main()->info("there are {} valid bonds after filtering", results.size());

    std::vector<chemical_entity::aminoacid const*> aminoacids;
    for (auto* res: _aminoacids)
        aminoacids.push_back(res);

    return {_pdb_name, params, aminoacids, results};
}
