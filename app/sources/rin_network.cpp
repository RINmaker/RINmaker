#include "rin_network.h"

#include "bonds.h"

#include "chemical_entity.h"

template<typename Bond>
void push_sort(std::list<const Bond*>& l, Bond const& b)
{
    static_assert(std::is_base_of<bond::base, Bond>::value, "template typename Bond must inherit from type bond");
    if (l.empty() || b < *l.front())
        l.push_front(&b);
    else
        l.push_back(&b);
}

void pairbond::push(bond::ss const& bond)
{ push_sort(_sss, bond); }

void pairbond::push(bond::vdw const& bond)
{ push_sort(_vdws, bond); }

void pairbond::push(bond::ionic const& bond)
{ push_sort(_ionics, bond); }

void pairbond::push(bond::hydrogen const& bond)
{ push_sort(_hydrogens, bond); }

void pairbond::push(bond::pication const& bond)
{ push_sort(_pications, bond); }

void pairbond::push(bond::pipistack const& bond)
{ push_sort(_pipistacks, bond); }

void pairbond::push(bond::generico const& bond)
{ push_sort(_generics, bond); }

pairbond::~pairbond()
{
    for (bond::base const* bond: _hydrogens)
        delete bond;

    for (bond::base const* bond: _sss)
        delete bond;

    for (bond::base const* bond: _vdws)
        delete bond;

    for (bond::base const* bond: _pications)
        delete bond;

    for (bond::base const* bond: _pipistacks)
        delete bond;

    for (bond::base const* bond: _ionics)
        delete bond;

    for (bond::base const* bond: _generics)
        delete bond;
}

network::~network()
{ for (auto const& kv: pairbonds_map) delete kv.second; }

pairbond& network::find(chemical_entity::aminoacid const& a, chemical_entity::aminoacid const& b)
{
    auto const key = prelude::concat_lexicographically(a.id(), b.id());
    if (pairbonds_map.find(key) == pairbonds_map.end())
    {
        auto* p = new pairbond();
        pairbonds_map.insert({key, p});
        return *p;
    }

    else
        return *pairbonds_map[key];
}

std::list<bond::base const*> network::get_one() const
{
    std::list<bond::base const*> tmp;
    for (auto const& kv: pairbonds_map)
    {
        auto* bond = kv.second->get_one();
        if (bond != nullptr)
            tmp.push_back(bond);
    }

    return tmp;
}

std::list<bond::base const*> network::get_all() const
{
    std::list<bond::base const*> tmp;
    for (auto const& kv: pairbonds_map)
    {
        auto lst = kv.second->get_all();
        tmp.splice(tmp.end(), lst);
    }

    return tmp;
}

std::list<bond::base const*> network::get_multiple() const
{
    std::list<bond::base const*> tmp;
    for (auto const& kv: pairbonds_map)
    {
        auto lst = kv.second->get_multiple();
        tmp.splice(tmp.end(), lst);
    }

    return tmp;
}


std::list<bond::base const*> pairbond::get_all() const
{
    std::list<bond::base const*> tmp;

    tmp.insert(tmp.end(), _hydrogens.begin(), _hydrogens.end());
    tmp.insert(tmp.end(), _sss.begin(), _sss.end());
    tmp.insert(tmp.end(), _vdws.begin(), _vdws.end());
    tmp.insert(tmp.end(), _pications.begin(), _pications.end());
    tmp.insert(tmp.end(), _pipistacks.begin(), _pipistacks.end());
    tmp.insert(tmp.end(), _ionics.begin(), _ionics.end());
    tmp.insert(tmp.end(), _generics.begin(), _generics.end());

    return tmp;
}

std::list<bond::base const*> pairbond::get_multiple() const
{
    std::list<bond::base const*> tmp;

    if (!_hydrogens.empty())
    {
        tmp.push_back(_hydrogens.front());
    }

    if (!_sss.empty())
    {
        tmp.push_back(_sss.front());
    }

    if (!_vdws.empty())
    {
        tmp.push_back(_vdws.front());
    }

    if (!_pications.empty())
    {
        tmp.push_back(_pications.front());
    }

    if (!_pipistacks.empty())
    {
        tmp.push_back(_pipistacks.front());
    }

    if (!_ionics.empty())
    {
        tmp.push_back(_ionics.front());
    }

    if (!_generics.empty())
    {
        tmp.push_back(_generics.front());
    }

    return tmp;
}

template<typename Bond>
bond::base const* __best(std::list<const Bond*> const& l, bond::base const* b)
{
    static_assert(std::is_base_of<bond::base, Bond>::value, "template typename Bond must inherit from type bond");
    if (!l.empty() && (b == nullptr || *b > *l.front()))
    {
        return l.front();
    }
    else
    {
        return b;
    }
}

bond::base const* pairbond::get_one() const
{
    bond::base const* best = nullptr;

    best = __best(_hydrogens, best);
    best = __best(_sss, best);
    best = __best(_vdws, best);
    best = __best(_ionics, best);
    best = __best(_pications, best);
    best = __best(_pipistacks, best);
    best = __best(_generics, best);

    return best;
}

std::list<bond::base const*> network::filter_hbond_realistic(std::list<bond::base const*> const& input) const
{
    std::set<bond::hydrogen const*> hydrogen_bonds_output;
    std::vector<bond::hydrogen const*> hydrogen_bonds_input;
    std::unordered_map<chemical_entity::atom const*, int> donors_bond_count;
    std::unordered_map<chemical_entity::atom const*, int> hydrogen_bond_count;
    std::unordered_map<chemical_entity::atom const*, int> acceptors_bond_count;

    //Get the bonds count of an atom
    auto get_bond_count = [](
            std::unordered_map<chemical_entity::atom const*, int>& container, chemical_entity::atom const* atom)->int
    {
        if (container.find(atom) == container.end())
            return 0;
        else
            return container[atom];
    };

    //Increase the bonds count of an atom
    auto inc_bond_count = [](
            std::unordered_map<chemical_entity::atom const*, int>& container, chemical_entity::atom const* atom)->void
    {
        if (container.find(atom) == container.end())
            container[atom] = 0;
        container[atom]++;
    };
    auto can_be_added = [&](bond::hydrogen const* bond)->bool
    {
        return (get_bond_count(donors_bond_count, bond->donor_ptr()) <
                bond->donor().how_many_hydrogen_can_donate() &&
                get_bond_count(hydrogen_bond_count, bond->hydrogen_ptr()) < 1 &&
                //An hydrogen can make only one bond
                get_bond_count(acceptors_bond_count, bond->acceptor_ptr()) <
                bond->acceptor().how_many_hydrogen_can_accept());
    };
    auto add_bond = [&](bond::hydrogen const* bond)->void
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
            hydrogen_bonds_input.push_back((bond::hydrogen const*) i);
    }

    //Order from smallest to largest energy
    sort(
            hydrogen_bonds_input.begin(), hydrogen_bonds_input.end(), [](bond::hydrogen const* a, bond::hydrogen const* b)
            { return a->get_energy() < b->get_energy(); });

    //Add as many hydrogen bonds as possible
    for (auto i: hydrogen_bonds_input)
    {
        if (can_be_added(i))
            add_bond(i);
    }

    //Let's build the output list
    std::list<bond::base const*> output;
    for (auto i: input)
    {
        //Insert i into the output if it is not an hydrogen or if it is in the filtered list
        if (i->get_type() != "hydrogen" ||
            hydrogen_bonds_output.find((bond::hydrogen const*) i) != hydrogen_bonds_output.end())
        {
            output.push_back(i);
        }
    }
    return output;
}