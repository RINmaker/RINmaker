#include "rin_graph.h"

#include <list>
#include <memory>
#include <queue>
#include <unordered_map>
#include <utility>

#include "chemical_entity.h"
#include "bonds.h"

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
    tmp_pimpl->source = bond.source_id();
    tmp_pimpl->target = bond.target_id();
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
    tmp_pimpl->source = bond.source_atom().res().id();
    tmp_pimpl->target = bond.target_atom().res().id();
    tmp_pimpl->distance = std::to_string(bond.get_length());
    tmp_pimpl->energy = std::to_string(bond.get_energy());
    tmp_pimpl->interaction = bond.get_interaction();
    tmp_pimpl->source_atom = bond.source_atom().name();
    tmp_pimpl->target_atom = bond.target_atom().name();
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
    tmp_pimpl->source = bond.source_positive().res().id();
    tmp_pimpl->target = bond.target_negative().res().id();
    tmp_pimpl->distance = std::to_string(bond.get_length());
    tmp_pimpl->energy = std::to_string(bond.get_energy());
    tmp_pimpl->interaction = bond.get_interaction();
    tmp_pimpl->source_atom = bond.source_positive().name();
    tmp_pimpl->target_atom = bond.target_negative().name();
    tmp_pimpl->positive = bond.source_positive().res().id();
    tmp_pimpl->angle = cfg::graphml::null;
    tmp_pimpl->donor = cfg::graphml::none;
    tmp_pimpl->cation = cfg::graphml::none;
    tmp_pimpl->orientation = cfg::graphml::none;
    pimpl = tmp_pimpl;
}

edge::edge(bond::hydrogen const& bond)
{
    auto tmp_pimpl = std::make_shared<impl>();
    tmp_pimpl->source = bond.source_atom().res().id();
    tmp_pimpl->target = bond.target_atom().res().id();
    tmp_pimpl->distance = std::to_string(bond.get_length());
    tmp_pimpl->energy = std::to_string(bond.get_energy());
    tmp_pimpl->interaction = bond.get_interaction();
    tmp_pimpl->source_atom = bond.source_atom().name();
    tmp_pimpl->target_atom = bond.target_atom().name();
    tmp_pimpl->angle = std::to_string(bond.get_angle());
    tmp_pimpl->donor = bond.donor().res().id();
    tmp_pimpl->cation = cfg::graphml::none;
    tmp_pimpl->positive = cfg::graphml::none;
    tmp_pimpl->orientation = cfg::graphml::none;
    pimpl = tmp_pimpl;
}

edge::edge(bond::pipistack const& bond)
{
    auto tmp_pimpl = std::make_shared<impl>();
    tmp_pimpl->source = bond.source_ring().res().id();
    tmp_pimpl->target = bond.target_ring().res().id();
    tmp_pimpl->distance = std::to_string(bond.get_length());
    tmp_pimpl->energy = std::to_string(bond.get_energy());
    tmp_pimpl->interaction = bond.get_interaction();
    tmp_pimpl->source_atom = bond.source_ring().name();
    tmp_pimpl->target_atom = bond.target_ring().name();
    tmp_pimpl->angle = std::to_string(bond.angle());
    tmp_pimpl->donor = cfg::graphml::none;
    tmp_pimpl->cation = cfg::graphml::none;
    tmp_pimpl->positive = cfg::graphml::none;
    tmp_pimpl->orientation = cfg::graphml::none;
    pimpl = tmp_pimpl;
}

edge::edge(bond::pication const& bond)
{
    auto tmp_pimpl = std::make_shared<impl>();
    tmp_pimpl->source = bond.source_ring().res().id();
    tmp_pimpl->target = bond.target_cation().res().id();
    tmp_pimpl->distance = std::to_string(bond.get_length());
    tmp_pimpl->energy = std::to_string(bond.get_energy());
    tmp_pimpl->interaction = bond.get_interaction();
    tmp_pimpl->source_atom = bond.source_ring().name();
    tmp_pimpl->target_atom = bond.target_cation().name();
    tmp_pimpl->cation = bond.target_cation().res().id();
    tmp_pimpl->angle = std::to_string(bond.angle());
    tmp_pimpl->donor = cfg::graphml::none;
    tmp_pimpl->positive = cfg::graphml::none;
    tmp_pimpl->orientation = cfg::graphml::none;
    pimpl = tmp_pimpl;
}

edge::edge(bond::generic_bond const& bond)
{
    auto tmp_pimpl = std::make_shared<impl>();
    tmp_pimpl->source = bond.source().id();
    tmp_pimpl->target = bond.target().id();
    tmp_pimpl->distance = std::to_string(bond.get_length());
    tmp_pimpl->energy = cfg::graphml::none;
    tmp_pimpl->interaction = bond.get_interaction();
    tmp_pimpl->source_atom = bond.source().name();
    tmp_pimpl->target_atom = bond.target().name();
    tmp_pimpl->angle = cfg::graphml::null;
    tmp_pimpl->donor = cfg::graphml::none;
    tmp_pimpl->cation = cfg::graphml::none;
    tmp_pimpl->positive = cfg::graphml::none;
    tmp_pimpl->orientation = cfg::graphml::none;
    pimpl = tmp_pimpl;
}

edge::~edge() = default;

string const& edge::source_id() const
{ return pimpl->source; }

string const& edge::target_id() const
{ return pimpl->target; }

string const& edge::distance() const
{ return pimpl->distance; }

string const& edge::energy() const
{ return pimpl->energy; }

string const& edge::interaction() const
{ return pimpl->interaction; }

string const& edge::source_atom() const
{ return pimpl->source_atom; }

string const& edge::target_atom() const
{ return pimpl->target_atom; }

string const& edge::angle() const
{ return pimpl->angle; }

string const& edge::donor() const
{ return pimpl->donor; }

string const& edge::cation() const
{ return pimpl->cation; }

string const& edge::positive() const
{ return pimpl->positive; }

string const& edge::orientation() const
{ return pimpl->orientation; }

void edge::append_to(xml_node& rin, bool with_metadata)
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
        parameters const& params,
        vector<aminoacid> const& aminoacids,
        vector<std::shared_ptr<bond::base const>> const& bonds)
{
    auto tmp_pimpl = std::make_shared<impl>(name, params);
    for (auto a: aminoacids)
    {
        auto n = (rin::node) a;
        tmp_pimpl->nodes.insert({n.get_id(), n});
    }

    // adjust nodes degree at edge insertion
    for (auto b: bonds)
    {
        auto edge = (rin::edge) *b;

        auto it = tmp_pimpl->nodes.find(edge.source_id());
        if (it != tmp_pimpl->nodes.end())
            it->second.inc_degree();

        it = tmp_pimpl->nodes.find(edge.target_id());
        if (it != tmp_pimpl->nodes.end())
            it->second.inc_degree();

        tmp_pimpl->edges.push_back(edge);
    }

    pimpl = tmp_pimpl;
}

graph::graph(graph const& other) : pimpl{new impl(*other.pimpl)}
{}

graph::~graph() = default;

string graph::name() const
{ return pimpl->name; }

vector<edge> graph::get_edges() const
{ return pimpl->edges; }

unordered_map<string, node> graph::get_nodes() const
{
    unordered_map<string, node> out;
    for (const auto& i: pimpl->nodes)
        out.insert_or_assign(i.first, i.second);

    return out;
}

void graph::write_to_file(fs::path const& out_path) const
{
    pugi::xml_document doc;

    // <graphml>
    pugi::xml_node graphml = doc.append_child("graphml");
    graphml.append_attribute("xmlns") = "http://graphml.graphdrawing.org/xmlns";
    graphml.append_attribute("xmlns:xsi") = "http://www.w3.org/2001/XMLSchema-instance";
    graphml.append_attribute("xsi:schemaLocation") =
            "http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd";

    // parameters summary
    auto comment = graphml.parent().insert_child_before(pugi::node_comment, graphml);
    comment.set_value(pimpl->params.pretty().c_str());

    // <graph>
    pugi::xml_node graph_node = graphml.append_child("graph");
    graph_node.append_attribute("id") = name().c_str();
    graph_node.append_attribute("edgedefault") = "undirected";

    // graphml requires all key attributes to be listed before the actual node/edges
    bool with_metadata = true;
    for (auto e: pimpl->edges)
    {
        e.append_to(graph_node, with_metadata);
        if (with_metadata) with_metadata = false;
    }

    with_metadata = true;
    for (auto kv: pimpl->nodes)
    {
        if (kv.second.degree() > 0)
        {
            kv.second.append_to(graph_node, with_metadata);
            if (with_metadata) with_metadata = false;
        }
    }

    doc.save_file(out_path.c_str());
}

node::node(chemical_entity::aminoacid const& res) : pimpl{new impl()}
{
    pimpl->id = res.id();
    pimpl->chain = res.chain_id();
    pimpl->seq = to_string(res.sequence_number());
    pimpl->name = res.name();
    pimpl->x = to_string(res[0]);
    pimpl->y = to_string(res[1]);
    pimpl->z = to_string(res[2]);
    pimpl->bfactor = res.ca() == nullptr ? "NULL" : to_string(res.ca()->temp_factor());
    pimpl->secondary = res.secondary_structure_id();
    pimpl->pdb_name = res.pdb_name();
    pimpl->degree = 0;
}

node::node(node const& other) : pimpl{new impl(*other.pimpl)}
{}

node& node::operator=(node const& rhs)
{
    if (this != &rhs)
        pimpl = std::make_unique<impl>(*rhs.pimpl);
    return *this;
}

node::~node() = default;

void node::inc_degree()
{ ++pimpl->degree; }

int node::degree() const
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
    add_data(pugi_node, "v_", "node", "Position", pimpl->seq, "double", with_metadata);
    add_data(pugi_node, "v_", "node", "Name", pimpl->name, "string", with_metadata);

    add_data(pugi_node, "v_", "node", "x", pimpl->x, "double", with_metadata);
    add_data(pugi_node, "v_", "node", "y", pimpl->y, "double", with_metadata);
    add_data(pugi_node, "v_", "node", "z", pimpl->z, "double", with_metadata);

    add_data(pugi_node, "v_", "node", "Bfactor_CA", pimpl->bfactor, "double", with_metadata);
    add_data(pugi_node, "v_", "node", "Secondary_Structure", pimpl->secondary, "string", with_metadata);

    add_data(pugi_node, "v_", "node", "PdbName", pimpl->pdb_name, "string", with_metadata);
}
