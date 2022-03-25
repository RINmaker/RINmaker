#include "rin_maker.h"

template<typename Record>
class secondary_structure_helper final {
private:
    std::unordered_map<std::string, std::map<interval<int>, Record, interval<int>::less>> _map;

public:
    void insert(Record const& record) {
        auto chain = _map.find(record.init_chain_id());
        if (chain == _map.end()) {
            std::map<interval<int>, Record, interval<int>::less> new_chain;
            new_chain.insert({record.range(), record});
            _map.insert({record.init_chain_id(), new_chain});
        } else {
            chain->second.insert({record.range(), record});
        }
    }

    [[nodiscard]]
    int size() const { return _map.size(); }

    // updates the secondary structure of res if res is present
    void update_if_contains(chemical_entity::aminoacid& res) const {
        auto chain = _map.find(res.chain_id());
        if (chain != _map.end()) {
            auto kv = chain->second.find(
                    interval<int>(res.sequence_number(), res.sequence_number()));
            if (kv != chain->second.end()) {
                res.make_secondary_structure(kv->second);
            }
        }
    }
};

rin_maker::base::base(std::filesystem::path const& pdb_path) {

    // open pdb
    std::ifstream pdb_file;
    pdb_file.open(pdb_path);

    // might throw
    if (!pdb_file.is_open())
        throw std::runtime_error("could not open " + pdb_path.string() + "\n");

    auto sheet_records = secondary_structure_helper<records::sheet_piece>();
    auto helix_records = secondary_structure_helper<records::helix>();

    std::vector<records::atom> tmp_atoms;

    log_manager::main()->info("parsing pdb...");

    // parsed line
    std::string line;
    while (getline(pdb_file, line)) {

        // first field is the type of the record
        std::string record_type = prelude::trim(line.substr(0, 6));

        if (record_type == "ATOM") {
            // atoms are grouped by aminoacid: we can simply fill a collection until we change residue
            records::atom record(line);
            if (!tmp_atoms.empty() && !record.same_res(tmp_atoms.back())) {
                _aminoacids.push_back(new chemical_entity::aminoacid(tmp_atoms));
                tmp_atoms.clear();
            }

            tmp_atoms.push_back(record);
        } else if (record_type == "HELIX") {
            helix_records.insert(records::helix(line));
        } else if (record_type == "SHEET") {
            sheet_records.insert(records::sheet_piece(line));
        } else if (record_type == "SSBOND") {
            _rin_network.new_bond<bonds::ss>(records::ss(line));
        }
    }

    if (!tmp_atoms.empty()) {
        _aminoacids.push_back(new chemical_entity::aminoacid(tmp_atoms));
    }

    log_manager::main()->info("done.");

    log_manager::main()->info(
            "finding appropriate secondary structure for amino acids...");

    if (sheet_records.size() != 0 || helix_records.size() != 0) {
        for (auto* res: _aminoacids) {
            res->make_secondary_structure();
            sheet_records.update_if_contains(*res);
            helix_records.update_if_contains(*res);
        }
    }

    log_manager::main()->info("done.");
}

template<typename BondFunc, typename Entity1, typename Entity2>
static void find_bonds(
        network& net, std::vector<Entity1 const*> const& vec, kdtree<Entity2, 3> const& tree, double distance) {
    static_assert(
            std::is_base_of<chemical_entity::component, Entity1>::value, "template typename Entity1 must inherit from type entity::component");
    static_assert(
            std::is_base_of<chemical_entity::component, Entity2>::value, "template typename Entity2 must inherit from type entity::component");
    static_assert(
            std::is_base_of<bondfunctors::base, BondFunc>::value, "template typename BondFunc must inherit from type bondfunctor::base");

    BondFunc if_test_insert(net);
    for (auto* e1: vec) {
        auto neighbors = tree.range_search(*e1, distance);
        for (auto* e2: neighbors)
            if_test_insert(*e1, *e2);
    }
}

rin_maker::all_bonds::all_bonds(std::filesystem::path const& pdb_path)
        : base(pdb_path) {

    log_manager::main()->info("retrieving components from aminoacids...");

    // used only to build the corresponding kdtrees
    std::vector<chemical_entity::atom const*> hdonors;
    std::vector<chemical_entity::ionic_group const*> positives;

    for (auto* res: _aminoacids) {

        for (auto const* a: res->atoms()) {
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
        if (ring != nullptr) {
            _ring_vector.push_back(ring);

            if (ring->is_a_pication_candidate())
                _pication_ring_vector.push_back(ring);
        }

        ring = res->secondary_ring();
        if (ring != nullptr) {
            _ring_vector.push_back(ring);

            if (ring->is_a_pication_candidate())
                _pication_ring_vector.push_back(ring);
        }
    }

    log_manager::main()->info("done.");

    log_manager::main()->info("building acceleration structures...");

    _hdonor_tree = kdtree<chemical_entity::atom, 3>(hdonors);
    _vdw_tree = kdtree<chemical_entity::atom, 3>(_vdw_vector);

    _ring_tree = kdtree<chemical_entity::ring, 3>(_ring_vector);
    _pication_ring_tree = kdtree<chemical_entity::ring, 3>(_pication_ring_vector);
    _positive_ion_tree = kdtree<chemical_entity::ionic_group, 3>(positives);

    log_manager::main()->info("done.");

    log_manager::main()->info("hydrogen acceptors: {}", _hacceptor_vector.size());
    log_manager::main()->info("hydrogen donors: {}", hdonors.size());
    log_manager::main()->info("vdw candidates: {}", _vdw_vector.size());
    log_manager::main()->info("cations: {}", _cation_vector.size());
    log_manager::main()->info("aromatic rings: {}", _ring_vector.size());
    log_manager::main()->info(
            "aromatic rings (valid for pication): {}", _ring_vector.size());

    find_bonds<bondfunctors::vdw>(
            _rin_network, _vdw_vector, _vdw_tree, parameters::get_distance_vdw() +
                                                  2 * cfg::params::max_vdw_radius);
    find_bonds<bondfunctors::ionic>(
            _rin_network, _negative_ion_vector, _positive_ion_tree, parameters::get_distance_ionic());
    find_bonds<bondfunctors::hydrogen>(
            _rin_network, _hacceptor_vector, _hdonor_tree, parameters::get_distance_h());
    find_bonds<bondfunctors::pication>(
            _rin_network, _cation_vector, _pication_ring_tree, parameters::get_distance_pication());
    find_bonds<bondfunctors::pipistack>(
            _rin_network, _ring_vector, _ring_tree, parameters::get_distance_pipistack());
}

rin_maker::backbone::backbone(std::filesystem::path const& pdb_path, carbon_getter getter)
        : base(pdb_path) {
    for (auto const* res: _aminoacids) {
        auto const* c = getter(*res);
        if (c != nullptr) {
            _carbon_vector.push_back(c);
        }
    }
    _carbon_tree = kdtree<chemical_entity::atom, 3>(_carbon_vector);

    find_bonds<bondfunctors::generico>(_rin_network, _carbon_vector, _carbon_tree, parameters::get_distance_generic());
}

rin::graph rin_maker::base::get_graph() const {
    rin::graph graph;

    for (auto* res: _aminoacids)
        graph.insert(res->to_node());

    std::list<bonds::base const*> results;
    switch (parameters::get_interaction_type()) {
        case parameters::interaction_type::MULTIPLE:
            results = _rin_network.get_multiple();
            break;

        case parameters::interaction_type::ALL:
            results = _rin_network.get_all();
            break;

        case parameters::interaction_type::ONE:
            results = _rin_network.get_one();
            break;
    }

    // TODO pensare ad un modo per mettere il filtro prima di get_multiple /
    // get_all ... (Può essere che tocchi creare una nuova net travasando solo i
    // risultati filtrati)
    if (parameters::get_hbond_realistic())
        results = _rin_network.filter_hbond_realistic(results);

    for (auto* bond: results)
        graph.push(bond->to_edge());

    log_manager::main()->info("there are {} bonds after filtering", results.size());

    return graph;
}

rin_maker::base const* rin_maker::build(parameters::policy policy, fs::path const& path)
{
    base const* rm = nullptr;
    switch (policy)
    {
        case parameters::policy::CLOSEST:
            rm = new all_bonds(path);
            break;
        case parameters::policy::CA:
            rm = new alpha_carbon(path);
            break;
        case parameters::policy::CB:
            rm = new beta_carbon(path);
            break;
    }

    return rm;
}