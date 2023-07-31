#include "rin_graph.h"

#include <list>
#include <memory>
#include <queue>
#include <unordered_map>

#include "ns_chemical_entity.h"
#include "ns_bond.h"

#include "config.h"

#include "private/impl_rin_graph.h"

namespace fs = std::filesystem;

using std::vector, std::string, std::queue, std::to_string, std::unordered_map, pugi::xml_node, std::to_string;
using chemical_entity::aminoacid;

using namespace rin;

void add_data(
    xml_node& parent,
    string const& prefix, string const& type,
    string const& key_name, string const& key_value, string const& key_type, bool with_metadata)
{
    xml_node data = parent.append_child("data");
    data.append_attribute("key") = (prefix + key_name).c_str();
    data.append_child(pugi::node_pcdata).set_value(key_value.c_str());

    if (with_metadata)
    {
        // add key info in (node -> graph -> graphml) before (node -> graph)
        xml_node key = parent.parent().parent().insert_child_before("key", parent.parent());
        key.append_attribute("id") = (prefix + key_name).c_str();
        key.append_attribute("for") = type.c_str();
        key.append_attribute("attr.name") = key_name.c_str();
        key.append_attribute("attr.type") = key_type.c_str();
    }
}


edge::edge(bond::ss const& bond)
{
    auto tmp_pimpl = std::make_shared<impl>();
    tmp_pimpl->source = bond.get_source_id();
    tmp_pimpl->target = bond.get_target_id();
    tmp_pimpl->distance = std::to_string(bond.get_length());
    tmp_pimpl->energy = std::to_string(bond.get_energy());
    tmp_pimpl->interaction = bond.get_interaction();
    tmp_pimpl->source_atom = "SG"; // TODO config
    tmp_pimpl->target_atom = "SG"; // TODO config
    tmp_pimpl->donor = cfg::graphml::none;
    tmp_pimpl->cation = cfg::graphml::none;
    tmp_pimpl->positive = cfg::graphml::none;
    tmp_pimpl->orientation = cfg::graphml::none;
    tmp_pimpl->angle = cfg::graphml::null;
    pimpl = tmp_pimpl;
}

edge::edge(bond::vdw const& bond)
{
    auto tmp_pimpl = std::make_shared<impl>();
    tmp_pimpl->source = bond.get_source_atom().get_residue().get_id();
    tmp_pimpl->target = bond.get_target_atom().get_residue().get_id();
    tmp_pimpl->distance = std::to_string(bond.get_length());
    tmp_pimpl->energy = std::to_string(bond.get_energy());
    tmp_pimpl->interaction = bond.get_interaction();
    tmp_pimpl->source_atom = bond.get_source_atom().get_name();
    tmp_pimpl->target_atom = bond.get_target_atom().get_name();
    tmp_pimpl->donor = cfg::graphml::none;
    tmp_pimpl->angle = cfg::graphml::null;
    tmp_pimpl->cation = cfg::graphml::none;
    tmp_pimpl->orientation = cfg::graphml::none;
    tmp_pimpl->positive = cfg::graphml::none;
    pimpl = tmp_pimpl;
}

edge::edge(bond::ionic const& bond)
{
    auto tmp_pimpl = std::make_shared<impl>();
    tmp_pimpl->source = bond.get_source_positive().get_residue().get_id();
    tmp_pimpl->target = bond.get_target_negative().get_residue().get_id();
    tmp_pimpl->distance = std::to_string(bond.get_length());
    tmp_pimpl->energy = std::to_string(bond.get_energy());
    tmp_pimpl->interaction = bond.get_interaction();
    tmp_pimpl->source_atom = bond.get_source_positive().get_name();
    tmp_pimpl->target_atom = bond.get_target_negative().get_name();
    tmp_pimpl->positive = bond.get_source_positive().get_residue().get_id();
    tmp_pimpl->angle = cfg::graphml::null;
    tmp_pimpl->donor = cfg::graphml::none;
    tmp_pimpl->cation = cfg::graphml::none;
    tmp_pimpl->orientation = cfg::graphml::none;
    pimpl = tmp_pimpl;
}

edge::edge(bond::hydrogen const& bond)
{
    auto tmp_pimpl = std::make_shared<impl>();
    tmp_pimpl->source = bond.get_source_atom().get_residue().get_id();
    tmp_pimpl->target = bond.get_target_atom().get_residue().get_id();
    tmp_pimpl->distance = std::to_string(bond.get_length());
    tmp_pimpl->energy = std::to_string(bond.get_energy());
    tmp_pimpl->interaction = bond.get_interaction();
    tmp_pimpl->source_atom = bond.get_source_atom().get_name();
    tmp_pimpl->target_atom = bond.get_target_atom().get_name();
    tmp_pimpl->angle = std::to_string(bond.get_angle());
    tmp_pimpl->donor = bond.get_donor().get_residue().get_id();
    tmp_pimpl->cation = cfg::graphml::none;
    tmp_pimpl->positive = cfg::graphml::none;
    tmp_pimpl->orientation = cfg::graphml::none;
    pimpl = tmp_pimpl;
}

edge::edge(bond::pipistack const& bond)
{
    auto tmp_pimpl = std::make_shared<impl>();
    tmp_pimpl->source = bond.get_source_ring().get_residue().get_id();
    tmp_pimpl->target = bond.get_target_ring().get_residue().get_id();
    tmp_pimpl->distance = std::to_string(bond.get_length());
    tmp_pimpl->energy = std::to_string(bond.get_energy());
    tmp_pimpl->interaction = bond.get_interaction();
    tmp_pimpl->source_atom = bond.get_source_ring().get_name();
    tmp_pimpl->target_atom = bond.get_target_ring().get_name();
    tmp_pimpl->angle = std::to_string(bond.get_angle());
    tmp_pimpl->donor = cfg::graphml::none;
    tmp_pimpl->cation = cfg::graphml::none;
    tmp_pimpl->positive = cfg::graphml::none;
    tmp_pimpl->orientation = cfg::graphml::none;
    pimpl = tmp_pimpl;
}

edge::edge(bond::pication const& bond)
{
    auto tmp_pimpl = std::make_shared<impl>();
    tmp_pimpl->source = bond.get_source_ring().get_residue().get_id();
    tmp_pimpl->target = bond.get_target_cation().get_residue().get_id();
    tmp_pimpl->distance = std::to_string(bond.get_length());
    tmp_pimpl->energy = std::to_string(bond.get_energy());
    tmp_pimpl->interaction = bond.get_interaction();
    tmp_pimpl->source_atom = bond.get_source_ring().get_name();
    tmp_pimpl->target_atom = bond.get_target_cation().get_name();
    tmp_pimpl->cation = bond.get_target_cation().get_residue().get_id();
    tmp_pimpl->angle = std::to_string(bond.get_angle());
    tmp_pimpl->donor = cfg::graphml::none;
    tmp_pimpl->positive = cfg::graphml::none;
    tmp_pimpl->orientation = cfg::graphml::none;
    pimpl = tmp_pimpl;
}

edge::edge(bond::generic_bond const& bond)
{
    auto tmp_pimpl = std::make_shared<impl>();
    tmp_pimpl->source = bond.get_source().get_id();
    tmp_pimpl->target = bond.get_target().get_id();
    tmp_pimpl->distance = std::to_string(bond.get_length());
    tmp_pimpl->energy = cfg::graphml::null;
    tmp_pimpl->interaction = bond.get_interaction();
    tmp_pimpl->source_atom = bond.get_source().get_name();
    tmp_pimpl->target_atom = bond.get_target().get_name();
    tmp_pimpl->angle = cfg::graphml::null;
    tmp_pimpl->donor = cfg::graphml::none;
    tmp_pimpl->cation = cfg::graphml::none;
    tmp_pimpl->positive = cfg::graphml::none;
    tmp_pimpl->orientation = cfg::graphml::none;
    pimpl = tmp_pimpl;
}

edge::~edge() = default;

string const& edge::get_source_id() const
{ return pimpl->source; }

string const& edge::get_target_id() const
{ return pimpl->target; }

string const& edge::get_distance() const
{ return pimpl->distance; }

string const& edge::get_energy() const
{ return pimpl->energy; }

string const& edge::get_interaction() const
{ return pimpl->interaction; }

string const& edge::get_source_atom() const
{ return pimpl->source_atom; }

string const& edge::get_target_atom() const
{ return pimpl->target_atom; }

string const& edge::get_angle() const
{ return pimpl->angle; }

string const& edge::get_donor() const
{ return pimpl->donor; }

string const& edge::get_cation() const
{ return pimpl->cation; }

string const& edge::get_positive() const
{ return pimpl->positive; }

string const& edge::get_orientation() const
{ return pimpl->orientation; }

void edge::append_to(xml_node& rin, bool with_metadata) const
{
    // the xml node representing a rin edge
    xml_node pugi_node = rin.append_child("edge");
    pugi_node.append_attribute("source") = pimpl->source.c_str();
    pugi_node.append_attribute("target") = pimpl->target.c_str();

    add_data(pugi_node, "e_", "edge", "NodeId1", pimpl->source, "string", with_metadata);
    add_data(pugi_node, "e_", "edge", "NodeId2", pimpl->target, "string", with_metadata);

    add_data(pugi_node, "e_", "edge", "Energy", pimpl->energy, "double", with_metadata);
    add_data(pugi_node, "e_", "edge", "Distance", pimpl->distance, "double", with_metadata);

    add_data(pugi_node, "e_", "edge", "Interaction", pimpl->interaction, "string", with_metadata);
    add_data(pugi_node, "e_", "edge", "Atom1", pimpl->source_atom, "string", with_metadata);
    add_data(pugi_node, "e_", "edge", "Atom2", pimpl->target_atom, "string", with_metadata);

    add_data(pugi_node, "e_", "edge", "Angle", pimpl->angle, "double", with_metadata);
    add_data(pugi_node, "e_", "edge", "Donor", pimpl->donor, "string", with_metadata);
    add_data(pugi_node, "e_", "edge", "Cation", pimpl->cation, "string", with_metadata);
    add_data(pugi_node, "e_", "edge", "Positive", pimpl->positive, "string", with_metadata);
    add_data(pugi_node, "e_", "edge", "Orientation", pimpl->orientation, "string", with_metadata);
}

graph::graph(
    string const& name,
    vector<aminoacid> const& aminoacids,
    vector<std::shared_ptr<bond::base const>> const& bonds)
{
    auto tmp_pimpl = std::make_shared<impl>(name);
    for (const auto& a: aminoacids)
    {
        auto n = (rin::node) a;
        tmp_pimpl->nodes.try_emplace(n.get_id(), n);
    }

    // adjust nodes degree at edge insertion
    for (const auto& b: bonds)
    {
        auto edge = (rin::edge) *b;

        auto it = tmp_pimpl->nodes.find(edge.get_source_id());
        if (it != tmp_pimpl->nodes.end())
            ++(it->second);

        it = tmp_pimpl->nodes.find(edge.get_target_id());
        if (it != tmp_pimpl->nodes.end())
            ++(it->second);

        tmp_pimpl->edges.push_back(edge);
    }

    pimpl = tmp_pimpl;
}

graph::graph(graph const& other) : pimpl{std::make_shared<impl>(*other.pimpl)}
{}

graph::~graph() = default;

string graph::get_name() const
{ return pimpl->name; }

vector<edge> const& graph::get_edges() const
{ return pimpl->edges; }

unordered_map<string, node> const& graph::get_nodes() const
{ return pimpl->nodes; }

void graph::write_to_file(fs::path const& out_path) const
{
    pugi::xml_document doc;

    // <graphml>
    pugi::xml_node graphml = doc.append_child("graphml");
    graphml.append_attribute("xmlns") = "http://graphml.graphdrawing.org/xmlns";
    graphml.append_attribute("xmlns:xsi") = "http://www.w3.org/2001/XMLSchema-instance";
    graphml.append_attribute("xsi:schemaLocation") =
            "http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd";

    // <graph>
    pugi::xml_node graph_node = graphml.append_child("graph");
    graph_node.append_attribute("id") = get_name().c_str();
    graph_node.append_attribute("edgedefault") = "undirected";

    // graphml requires all key attributes to be listed before the actual node/edges
    bool first_time = true;
    for (const auto& e : pimpl->edges)
    {
        if (first_time)
        {
            e.append_to(graph_node, true);
            first_time = false;
        }
        else
            e.append_to(graph_node);
    }

    first_time = true;
    for (auto const& [_, node] : pimpl->nodes)
    {
        if (node.get_degree() > 0)
        {
            if (first_time)
            {
                node.append_to(graph_node, true);
                first_time = false;
            }
            else
                node.append_to(graph_node);
        }
    }

    doc.save_file(out_path.c_str());
}

node::node(chemical_entity::aminoacid const& res) : pimpl{std::make_unique<impl>()}
{
    pimpl->id = res.get_id();
    pimpl->chain = res.get_chain_id();
    pimpl->sequence_number = to_string(res.get_sequence_number());
    pimpl->name = res.get_name();
    pimpl->x = to_string(res.get_position()[0]);
    pimpl->y = to_string(res.get_position()[1]);
    pimpl->z = to_string(res.get_position()[2]);
    pimpl->bfactor_ca = res.get_alpha_carbon().has_value() ? to_string(res.get_alpha_carbon().value().get_temp_factor()) : cfg::graphml::null;
    pimpl->secondary_structure = res.get_secondary_structure_id();
    pimpl->pdb_name = res.get_protein_name();
    pimpl->degree = 0;
}

node::node(node const& other) : pimpl{std::make_unique<impl>(*other.pimpl)}
{}

node& node::operator=(node const& rhs)
{
    if (this != &rhs)
        pimpl = std::make_unique<impl>(*rhs.pimpl);
    return *this;
}

node::~node() = default;

node& node::operator++()
{ ++pimpl->degree; return *this; }

int node::get_degree() const
{ return pimpl->degree; }

[[nodiscard]]
string const& node::get_id() const
{ return pimpl->id; }

void node::append_to(xml_node& graph, bool with_metadata) const
{
    xml_node pugi_node;

    pugi_node = graph.prepend_child("node");
    pugi_node.append_attribute("id") = pimpl->id.c_str();

    add_data(pugi_node, "v_", "node", "Degree", to_string(pimpl->degree), "double", with_metadata);
    add_data(pugi_node, "v_", "node", "NodeId", pimpl->id, "string", with_metadata);

    add_data(pugi_node, "v_", "node", "Residue", pimpl->id, "string", with_metadata);
    add_data(pugi_node, "v_", "node", "Chain", pimpl->chain, "string", with_metadata);
    add_data(pugi_node, "v_", "node", "Position", pimpl->sequence_number, "double", with_metadata);
    add_data(pugi_node, "v_", "node", "Name", pimpl->name, "string", with_metadata);

    add_data(pugi_node, "v_", "node", "x", pimpl->x, "double", with_metadata);
    add_data(pugi_node, "v_", "node", "y", pimpl->y, "double", with_metadata);
    add_data(pugi_node, "v_", "node", "z", pimpl->z, "double", with_metadata);

    add_data(pugi_node, "v_", "node", "Bfactor_CA", pimpl->bfactor_ca, "double", with_metadata);
    add_data(pugi_node, "v_", "node", "Secondary_Structure", pimpl->secondary_structure, "string", with_metadata);

    add_data(pugi_node, "v_", "node", "PdbName", pimpl->pdb_name, "string", with_metadata);
}
