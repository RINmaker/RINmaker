#include "rin_graph.h"

#include <filesystem>
#include <string>

#include "noncovalent_bonds.h"
#include "chemical_entity.h"
#include "config.h"

using namespace rin;
using chemical_entity::aminoacid;

graph::graph(rin::parameters const& params, vector<aminoacid const*> const& aminoacids, list<bonds::base const*> const& bonds) : _params(params)
{
    for (auto a: aminoacids)
    {
        auto n = a->to_node();
        _nodes.insert({n.get_id(), n});
    }

    // adjust _nodes degree at edge insertion
    for (auto b: bonds)
    {
        auto edge = (rin::edge) *b;

        auto it = _nodes.find(edge.source_id());
        if (it != _nodes.end())
            ++(it->second.degree());

        it = _nodes.find(edge.target_id());
        if (it != _nodes.end())
            ++(it->second.degree());

        _edges.push_back(edge);
    }
}

void graph::write_to_file(std::filesystem::path const& out_path)
{
    pugi::xml_document doc;

    // <graphml>
    pugi::xml_node graphml = doc.append_child("graphml");
    graphml.append_attribute("xmlns") = "http://graphml.graphdrawing.org/xmlns";
    graphml.append_attribute("xmlns:xsi") = "http://www.w3.org/2001/XMLSchema-instance";
    graphml.append_attribute("xsi:schemaLocation") = "http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd";

    // parameters summary
    auto comment = graphml.parent().insert_child_before(pugi::node_comment, graphml);
    comment.set_value(_params.pretty().c_str());

    // <graph>
    pugi::xml_node graph_node = graphml.append_child("graph");
    graph_node.append_attribute("id") = "G";
    graph_node.append_attribute("edgedefault") = "undirected";

    // graphml requires all key attributes to be listed before the actual node/edges
    bool with_metadata = true;
    for (auto e: _edges)
    {
        e.append_to(graph_node, with_metadata);
        if (with_metadata) with_metadata = false;
    }

    with_metadata = true;
    for (auto kv: _nodes)
    {
        if (kv.second.degree() > 0)
        {
            kv.second.append_to(graph_node, with_metadata);
            if (with_metadata) with_metadata = false;
        }
    }

    doc.save_file(out_path.c_str());
}

edge::edge(bonds::ss const& bond) :
        _source(bond.source_id()),
        _target(bond.target_id()),
        _distance(std::to_string(bond.get_length())),
        _energy(std::to_string(bond.get_energy())),
        _interaction(bond.get_interaction()),
        _source_atom("SG"), // TODO va in config
        _target_atom("SG"), // TODO va in config
        _donor(cfg::graphml::none),
        _cation(cfg::graphml::none),
        _positive(cfg::graphml::none),
        _orientation(cfg::graphml::none),
        _angle(cfg::graphml::null)
{}

edge::edge(bonds::vdw const& bond) :
        _source(bond.source_atom().res().id()),
        _target(bond.target_atom().res().id()),
        _distance(std::to_string(bond.get_length())),
        _energy(std::to_string(bond.get_energy())),
        _interaction(bond.get_interaction()),
        _source_atom(bond.source_atom().name()),
        _target_atom(bond.target_atom().name()),
        _donor(cfg::graphml::none),
        _angle(cfg::graphml::null),
        _cation(cfg::graphml::none),
        _orientation(cfg::graphml::none),
        _positive(cfg::graphml::none)
{}

edge::edge(bonds::ionic const& bond) :
        _source(bond.positive().res().id()),
        _target(bond.negative().res().id()),
        _distance(std::to_string(bond.get_length())),
        _energy(std::to_string(bond.get_energy())),
        _interaction(bond.get_interaction()),
        _source_atom(bond.positive().name()),
        _target_atom(bond.negative().name()),
        _positive(bond.positive().res().id()),
        _angle(cfg::graphml::null),
        _donor(cfg::graphml::none),
        _cation(cfg::graphml::none),
        _orientation(cfg::graphml::none)
{}

edge::edge(bonds::hydrogen const& bond)
        :
        _source(bond.acceptor().res().id()),
        _target(bond.donor().res().id()),
        _distance(std::to_string(bond.get_length())),
        _energy(std::to_string(bond.get_energy())),
        _interaction(bond.get_interaction()),
        _source_atom(bond.acceptor().name()),
        _target_atom(bond.donor().name()),
        _angle(std::to_string(bond.get_angle())),
        _donor(bond.donor().res().id()),
        _cation(cfg::graphml::none),
        _positive(cfg::graphml::none),
        _orientation(cfg::graphml::none)
{}

edge::edge(bonds::pipistack const& bond) :
        _source(bond.source_ring().res().id()),
        _target(bond.target_ring().res().id()),
        _distance(std::to_string(bond.get_length())),
        _energy(std::to_string(bond.get_energy())),
        _interaction(bond.get_interaction()),
        _source_atom(bond.source_ring().name()),
        _target_atom(bond.target_ring().name()),
        _angle(std::to_string(bond.angle())),
        _donor(cfg::graphml::none),
        _cation(cfg::graphml::none),
        _positive(cfg::graphml::none),
        _orientation(cfg::graphml::none)
{}

edge::edge(bonds::pication const& bond) :
        _source(bond.ring().res().id()),
        _target(bond.cation().res().id()),
        _distance(std::to_string(bond.get_length())),
        _energy(std::to_string(bond.get_energy())),
        _interaction(bond.get_interaction()),
        _source_atom(bond.ring().name()),
        _target_atom(bond.cation().name()),
        _cation(bond.cation().res().id()),
        _angle(std::to_string(bond.angle())),
        _donor(cfg::graphml::none),
        _positive(cfg::graphml::none),
        _orientation(cfg::graphml::none)
{}

edge::edge(bonds::generico const& bond)
        :
        _source(bond.source().id()),
        _target(bond.target().id()),
        _distance(std::to_string(bond.get_length())),
        _energy(cfg::graphml::none),
        _interaction(bond.get_interaction()),
        _source_atom(bond.source().name()),
        _target_atom(bond.target().name()),
        _angle(cfg::graphml::null),
        _donor(cfg::graphml::none),
        _cation(cfg::graphml::none),
        _positive(cfg::graphml::none),
        _orientation(cfg::graphml::none)
{}

void add_data(
        pugi::xml_node& node, string const& prefix, string const& type, string const& key_name, string const& key_value, string const& key_type, bool with_metadata)
{
    pugi::xml_node data = node.append_child("data");
    data.append_attribute("key") = (prefix + key_name).c_str();
    data.append_child(pugi::node_pcdata).set_value(key_value.c_str());

    if (with_metadata)
    {
        // add key info in (node -> graph -> graphml) before (node -> graph)
        pugi::xml_node key = node.parent().parent().insert_child_before("key", node.parent());
        key.append_attribute("id") = (prefix + key_name).c_str();
        key.append_attribute("for") = type.c_str();
        key.append_attribute("attr.name") = key_name.c_str();
        key.append_attribute("attr.type") = key_type.c_str();
    }
}

void edge::append_to(pugi::xml_node& rin, bool metadata)
{
    // the xml node representing a rin edge
    pugi::xml_node edge = rin.append_child("edge");
    edge.append_attribute("source") = _source.c_str();
    edge.append_attribute("target") = _target.c_str();

    add_data(edge, "e_", "edge", "NodeId1", _source, "string", metadata);
    add_data(edge, "e_", "edge", "NodeId2", _target, "string", metadata);

    add_data(edge, "e_", "edge", "Energy", _energy, "double", metadata);
    add_data(edge, "e_", "edge", "Distance", _distance, "double", metadata);

    add_data(edge, "e_", "edge", "Interaction", _interaction, "string", metadata);
    add_data(edge, "e_", "edge", "Atom1", _source_atom, "string", metadata);
    add_data(edge, "e_", "edge", "Atom2", _target_atom, "string", metadata);

    add_data(edge, "e_", "edge", "Angle", _angle, "double", metadata);
    add_data(edge, "e_", "edge", "Donor", _donor, "string", metadata);
    add_data(edge, "e_", "edge", "Cation", _cation, "string", metadata);
    add_data(edge, "e_", "edge", "Positive", _positive, "string", metadata);
    add_data(edge, "e_", "edge", "Orientation", _orientation, "string", metadata);
}

node::node(chemical_entity::aminoacid const& res) :
        _id(res.id()),
        _chain(res.chain_id()),
        _seq(std::to_string(res.sequence_number())),
        _name(res.name()),
        _x(std::to_string(res[0])),
        _y(std::to_string(res[1])),
        _z(std::to_string(res[2])),
        _bfactor(res.ca() == nullptr ? "NULL" : std::to_string(res.ca()->temp_factor())),
        _secondary(res.secondary_structure_id()),
        _pdb_name(res.pdb_name()), // TODO
        _degree(0)
{}


void node::append_to(pugi::xml_node& graph, bool with_metadata) const
{
    pugi::xml_node node;

    node = graph.prepend_child("node");
    node.append_attribute("id") = _id.c_str();

    add_data(node, "v_", "node", "Degree", std::to_string(_degree), "double", with_metadata);
    add_data(node, "v_", "node", "NodeId", _id, "string", with_metadata);

    add_data(node, "v_", "node", "Residue", _id, "string", with_metadata);
    add_data(node, "v_", "node", "Chain", _chain, "string", with_metadata);
    add_data(node, "v_", "node", "Position", _seq, "double", with_metadata);
    add_data(node, "v_", "node", "Name", _name, "string", with_metadata);

    add_data(node, "v_", "node", "x", _x, "double", with_metadata);
    add_data(node, "v_", "node", "y", _y, "double", with_metadata);
    add_data(node, "v_", "node", "z", _z, "double", with_metadata);

    add_data(node, "v_", "node", "Bfactor_CA", _bfactor, "double", with_metadata);
    add_data(node, "v_", "node", "Secondary_Structure", _secondary, "string", with_metadata);

    add_data(node, "v_", "node", "PdbName", _pdb_name, "string", with_metadata);
}
