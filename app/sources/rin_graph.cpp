#include "rin_graph.h"

#include <filesystem>
#include <string>

#include "noncovalent_bonds.h"
#include "chemical_entity.h"
#include "config.h"

using namespace rin;

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

// add <data/> to a pugixml node
//
void __add_edge_data(pugi::xml_node& node, string key_name, string key_value, string key_type, bool metadata)
{
    // add data to node
    pugi::xml_node data = node.append_child("data");
    data.append_attribute("key") = ("e_" + key_name).c_str();
    data.append_child(pugi::node_pcdata).set_value(key_value.c_str());

    if (metadata)
    {
        // add key info in (node -> graph -> graphml) before (node -> graph)
        pugi::xml_node key = node.parent().parent().insert_child_before("key", node.parent());
        key.append_attribute("id") = ("e_" + key_name).c_str();
        key.append_attribute("for") = "edge";
        key.append_attribute("attr.name") = key_name.c_str();
        key.append_attribute("attr.type") = key_type.c_str();
    }
}

void edge::append_to(pugi::xml_node& rin, bool metadata)
{
    // the xml node representing a rin edge
    pugi::xml_node e = rin.append_child("edge");
    e.append_attribute("source") = _source.c_str();
    e.append_attribute("target") = _target.c_str();

    __add_edge_data(e, "NodeId1", _source, "string", metadata);
    __add_edge_data(e, "NodeId2", _target, "string", metadata);

    __add_edge_data(e, "Energy", _energy, "double", metadata);
    __add_edge_data(e, "Distance", _distance, "double", metadata);

    __add_edge_data(e, "Interaction", _interaction, "string", metadata);
    __add_edge_data(e, "Atom1", _source_atom, "string", metadata);
    __add_edge_data(e, "Atom2", _target_atom, "string", metadata);

    __add_edge_data(e, "Angle", _angle, "double", metadata);
    __add_edge_data(e, "Donor", _donor, "string", metadata);
    __add_edge_data(e, "Cation", _cation, "string", metadata);
    __add_edge_data(e, "Positive", _positive, "string", metadata);
    __add_edge_data(e, "Orientation", _orientation, "string", metadata);
}

edge graph::pop_edge()
{
    edge e = edges.front();

    auto it = nodes.find(e.source_id());
    if (it != nodes.end())
    {
        ++(it->second.degree());
    }

    it = nodes.find(e.target_id());
    if (it != nodes.end())
    {
        ++(it->second.degree());
    }

    edges.pop();
    return e;
}

void graph::consume_to_xml(rin::parameters const& params, std::filesystem::path const& out_path)
{
    pugi::xml_document doc;

    // <graphml>
    pugi::xml_node graphml = doc.append_child("graphml");
    graphml.append_attribute("xmlns") = "http://graphml.graphdrawing.org/xmlns";
    graphml.append_attribute("xmlns:xsi") = "http://www.w3.org/2001/XMLSchema-instance";
    graphml.append_attribute("xsi:schemaLocation") = "http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd";

    // parameters summary
    auto comment = graphml.parent().insert_child_before(pugi::node_comment, graphml);
    comment.set_value(params.pretty().c_str());

    // <graph>
    pugi::xml_node graph_node = graphml.append_child("graph");
    graph_node.append_attribute("id") = "G";
    graph_node.append_attribute("edgedefault") = "undirected";

    // graphml requires all key attributes to be listed before the actual node/edges
    bool metadata = true;
    while (!edges.empty())
    {
        edge e = pop_edge();
        e.append_to(graph_node, metadata);
        if (metadata) metadata = false;
    }

    metadata = true;
    for (auto kv: nodes)
    {
        if (kv.second.degree() > 0)
        {
            kv.second.append_to(graph_node, metadata);
            if (metadata)
            {
                metadata = false;
            }
        }
    }

    // if (!out_path.has_parent_path())
    // std::filesystem::create_directory(out_path.parent_path());

    doc.save_file(out_path.c_str());
}

node::node(chemical_entity::aminoacid const& res)
        :
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

// add <data/> to a pugixml node
//
void __add_node_data(pugi::xml_node& node, std::string key_name, std::string key_value, std::string key_type, bool metadata)
{
    // add data to node
    pugi::xml_node data = node.append_child("data");
    data.append_attribute("key") = ("v_" + key_name).c_str();
    data.append_child(pugi::node_pcdata).set_value(key_value.c_str());

    if (metadata)
    {
        // add key info in (node -> graph -> graphml) before (node -> graph)
        pugi::xml_node key = node.parent().parent().insert_child_before("key", node.parent());
        key.append_attribute("id") = ("v_" + key_name).c_str();
        key.append_attribute("for") = "edge";
        key.append_attribute("attr.name") = key_name.c_str();
        key.append_attribute("attr.type") = key_type.c_str();
    }
}

void node::append_to(pugi::xml_node& graph, bool metadata) const
{
    pugi::xml_node node;

    node = graph.prepend_child("node");
    node.append_attribute("id") = _id.c_str();
    __add_node_data(node, "Degree", std::to_string(_degree), "double", metadata);
    __add_node_data(node, "NodeId", _id, "string", metadata);

    __add_node_data(node, "Residue", _id, "string", metadata);
    __add_node_data(node, "Chain", _chain, "string", metadata);
    __add_node_data(node, "Position", _seq, "double", metadata);
    __add_node_data(node, "Name", _name, "string", metadata);

    __add_node_data(node, "x", _x, "double", metadata);
    __add_node_data(node, "y", _y, "double", metadata);
    __add_node_data(node, "z", _z, "double", metadata);

    __add_node_data(node, "Bfactor_CA", _bfactor, "double", metadata);
    __add_node_data(node, "Secondary_Structure", _secondary, "string", metadata);

    __add_node_data(node, "PdbName", _pdb_name, "string", metadata);
}


std::vector<edge> graph::get_edges()
{
    std::queue<edge> q(edges);
    std::vector<edge> out;
    for (int i = 0; !q.empty(); i++)
    {
        out.push_back(q.front());
        q.pop();
    }
    return out;
}

std::unordered_map<std::string, node> graph::get_nodes()
{
    std::unordered_map<std::string, node> out;
    for (const auto& i: nodes)
    {
        out.insert_or_assign(i.first, i.second);
    }
    return out;
}