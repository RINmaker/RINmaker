#include "pdb_data.h"

#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <exception>

#include "bonds.h"
#include "records.h"
#include "kdtree.h"
#include "config.h"
#include "prelude.h"
#include "interval.h"
#include "bondfunctors.h"

using namespace std;

template<typename BondFunc, typename Entity1, typename Entity2>
void __find_bonds(network& net, std::vector<Entity1 const*> const& vec, kdtree<Entity2, 3> const& tree, double distance)
{
    static_assert(std::is_base_of<entities::component, Entity1>::value, "template typename Entity1 must inherit from type component");
    static_assert(std::is_base_of<entities::component, Entity2>::value, "template typename Entity2 must inherit from type component");
    static_assert(std::is_base_of<bondfunctors::base, BondFunc>::value, "template typename BondFunc must inherit from type bondfunctors::base");

    BondFunc if_test_insert(net);
    for (auto* e1: vec)
    {
        auto neighbors = tree.range_search(*e1, distance);
        for (auto* e2: neighbors)
            if_test_insert(*e1, *e2);
    }
}

pdb_data::pdb_data(bool exportOutput) : logger(log_manager::main())
{
    try
    {
        // might throw
        log_manager::main()->info(">> begin raw pdb data <<");
        parse_file();
        log_manager::main()->info(">> end raw pdb data <<");

        // compute bonds
        log_manager::main()->info(">> begin computed bond data <<");
        switch (parameters::get_net_policy())
        {
            case parameters::policy::CLOSEST:
                // è leggermente diversa per come loro vogliono specificare questa distanza... non vogliono specificare
                // quella tra i centri (come in tutti gli altri) ma quella tra le superfici delle sfere, che è comunque
                // relata a quella tra i centri in modo semplice.
                //
                __find_bonds<bondfunctors::vdw>(net, _vdws, _vdws_tree, parameters::get_distance_vdw() + 2 * cfg::params::max_vdw_radius);

                __find_bonds<bondfunctors::ionic>(net, _negatives, _positives_tree, parameters::get_distance_ionic());
                __find_bonds<bondfunctors::hydrogen>(net, _hacceptors, _hdonors_tree, parameters::get_distance_h());
                __find_bonds<bondfunctors::pication>(net, _cations, _pication_rings_tree, parameters::get_distance_pication());
                __find_bonds<bondfunctors::pipistack>(net, _rings, _rings_tree, parameters::get_distance_pipistack());
                break;

            case parameters::policy::CA:
                __find_bonds<bondfunctors::generico>(net, _cas, _cas_tree, parameters::get_distance_generic());
                break;

            case parameters::policy::CB:
                __find_bonds<bondfunctors::generico>(net, _cbs, _cbs_tree, parameters::get_distance_generic());
                break;
        }

        // scrivi in grafo e to xml
        for (auto* res: _aminoacids)
            g.insert(res->to_node());

        list<bonds::base const*> results;
        switch (parameters::get_interaction_type())
        {
            case parameters::interaction_type::MULTIPLE:
                results = net.get_multiple();
                break;

            case parameters::interaction_type::ALL:
                results = net.get_all();
                break;

            case parameters::interaction_type::ONE:
                results = net.get_one();
                break;
        }

        //TODO pensare ad un modo per mettere il filtro prima di get_multiple / get_all ... (Può essere che tocchi creare una nuova net travasando solo i risultati filtrati)
        if(parameters::get_hbond_realistic())
            results = net.filter_hbond_realistic(results);

        for (auto* bond: results)
            g.push(bond->to_edge());

        log_manager::main()->info("there are {} bonds after filtering", results.size());
        log_manager::main()->info(">> end computed bond data <<");

        // write to output
        if (exportOutput)
        {
            g.consume_to_xml();
        }
    }

    catch (runtime_error& e)
    {
        logger->error(e.what());
        throw;
    }
}

void pdb_data::parse_file()
{
    // aprire il pdb

    ifstream pdb_file;
    pdb_file.open(parameters::get_pdb_path());

    if (!pdb_file.is_open())
    {
        throw runtime_error("could not open " + parameters::get_pdb_path().string() + "\n");
    }

    // parsare il file pdb

    // organizzare gli sheet in una struttura sortata, divisi per catene
    unordered_map<string, map<interval<int>, records::sheet_piece, interval<int>::less>> catene_sheets;

    // organizzare le helix in una struttura sortata, divise per catene
    unordered_map<string, map<interval<int>, records::helix, interval<int>::less>> catene_helices;

    vector<records::atom> tmp_atoms;

    // parsed line
    string line;
    while (getline(pdb_file, line))
    {
        // first field is the type of the record
        string intestazione = prelude::trim(line.substr(0, 6));


        if (intestazione == "ATOM")
        {
            // atoms are grouped by aminoacid: we can simply fill a collection until we change residue

            records::atom record(line);
            if (!tmp_atoms.empty() && !record.same_res(tmp_atoms.back()))
            {
                _aminoacids.push_back(new entities::aminoacid(tmp_atoms));
                tmp_atoms.clear();
            }
            tmp_atoms.push_back(record);
        }
        else if (intestazione == "HELIX")
        {
            // inserisco la helix direttamente nella struttura divisa per catene e sortata sugli intervalli di amminoacidi

            records::helix record(line);
            auto catena = catene_helices.find(record.init_chain_id());
            if (catena == catene_helices.end())
            {
                map<interval<int>, records::helix, interval<int>::less> nuova_catena;
                nuova_catena.insert({record.range(), record});
                catene_helices.insert({record.init_chain_id(), nuova_catena});
            }

            else
            {
                catena->second.insert({record.range(), record});
            }
        }

            // inserisco lo sheet direttamente nella struttura divisa per catene e sortata sugli intervalli di amminoacidi
        else if (intestazione == "SHEET")
        {
            records::sheet_piece record(line);

            auto catena = catene_sheets.find(record.init_chain_id());
            if (catena == catene_sheets.end())
            {
                map<interval<int>, records::sheet_piece, interval<int>::less> nuova_catena;
                nuova_catena.insert({record.range(), record});
                catene_sheets.insert({record.init_chain_id(), nuova_catena});
            }

            else
            {
                catena->second.insert({record.range(), record});
            }
        }

            // per l'ssbond invece lo inserisco diretto nella rete
        else if (intestazione == "SSBOND")
        {
            net.new_bond<bonds::ss>(records::ss(line));
        }
    }

    // ultimo controllo sui record di amminoacido
    if (!tmp_atoms.empty())
    {
        _aminoacids.push_back(new entities::aminoacid(tmp_atoms));
    }

    // aggiustare gli amminoacidi e recuperare i vari componenti
    std::vector<entities::atom const*> hdonors;
    std::vector<entities::ionic_group const*> positives;

    for (auto* res: _aminoacids)
    {
        auto& resr = *res;

        // recuperare i singoli atomi sse candidati a legami
        for (auto const* a: resr.atoms())
        {
            auto& ar = *a;
            if (ar.is_a_hydrogen_donor())
            {
                hdonors.push_back(a);
            }

            if (ar.is_a_hydrogen_acceptor())
            {
                _hacceptors.push_back(a);
            }

            if (ar.is_a_vdw_candidate())
            {
                _vdws.push_back(a);
            }

            if (ar.is_a_cation())
            {
                _cations.push_back(a);
            }
        }

        // recuperare i gruppi ionici
        entities::ionic_group const* group = resr.positive_ionic_group();
        if (group != nullptr)
        {
            positives.push_back(group);
        }

        group = resr.negative_ionic_group();
        if (group != nullptr)
        {
            _negatives.push_back(group);
        }

        // recuperare i ring
        entities::ring const* ring = resr.primary_ring();
        if (ring != nullptr)
        {
            _rings.push_back(ring);
            if (ring->is_a_pication_candidate())
            {
                _pication_rings.push_back(ring);
            }
        }

        ring = resr.secondary_ring();
        if (ring != nullptr)
        {
            _rings.push_back(ring);
            if (ring->is_a_pication_candidate())
            {
                _pication_rings.push_back(ring);
            }
        }

        // recuperare alpha e beta carbon
        auto const* c = resr.ca();
        if (c != nullptr)
        {
            _cas.push_back(c);
        }

        c = resr.cb();
        if (c != nullptr)
        {
            _cbs.push_back(c);
        }

        // stabilire se l'amminoacido è in una helix, sheet (oppure nulla)
        // in questo punto sappiamo se il file contiene informazioni oppure no!
        //
        if (catene_sheets.size() != 0 || catene_helices.size() != 0)
        {
            // setto una struttura di base intanto
            resr.make_secondary_structure();

            // è in uno sheet?
            auto catena_sheet = catene_sheets.find(res->chain_id());
            if (catena_sheet != catene_sheets.end())
            {
                auto kv = catena_sheet->second.find(interval<int>(res->sequence_number(), res->sequence_number()));
                if (kv != catena_sheet->second.end())
                {
                    resr.make_secondary_structure(kv->second);
                }
            }

            // oppure è in una helix?
            auto catena_helix = catene_helices.find(res->chain_id());
            if (catena_helix != catene_helices.end())
            {
                auto kv = catena_helix->second.find(interval<int>(res->sequence_number(), res->sequence_number()));
                if (kv != catena_helix->second.end())
                {
                    resr.make_secondary_structure(kv->second);
                }
            }
        }
    }

    log_manager::main()->info("alpha carbons: {}", _cas.size());
    log_manager::main()->info("beta carbons: {}", _cbs.size());

    log_manager::main()->info("hydrogen acceptors: {}", _hacceptors.size());
    log_manager::main()->info("hydrogen donors: {}", hdonors.size());
    log_manager::main()->info("vdw candidates: {}", _vdws.size());
    log_manager::main()->info("cations: {}", _cations.size());
    log_manager::main()->info("aromatic rings: {}", _rings.size());
    log_manager::main()->info("aromatic rings (valid for pication): {}", _rings.size());

    // costruire gli alberi dei donatori di idrogeno e dei possibili partecipanti in van der waals
    _hdonors_tree = kdtree<entities::atom, 3>(hdonors);
    _vdws_tree = kdtree<entities::atom, 3>(_vdws);

    // costruire gli alberi dei ring (pipistack e pication) e dei gruppi ionici (ionic bond)
    _rings_tree = kdtree<entities::ring, 3>(_rings);
    _pication_rings_tree = kdtree<entities::ring, 3>(_pication_rings);
    _positives_tree = kdtree<entities::ionic_group, 3>(positives);

    // costruire gli alberi dei carboni alpha e beta
    _cas_tree = kdtree<entities::atom, 3>(_cas);
    _cbs_tree = kdtree<entities::atom, 3>(_cbs);
}