#include "rin_network.h"

#include "bonds.h"

#include "chemical_entity.h"

template<typename Bond>
void push_sort(std::list<std::shared_ptr<Bond const>>& l, std::shared_ptr<Bond const> b)
{
    static_assert(std::is_base_of<bond::base, Bond>::value, "template typename Bond must inherit from type bond::base");
    if (l.empty() || *b < *l.front())
        l.push_front(b);
    else
        l.push_back(b);
}

using std::shared_ptr;

void pairbond::push(shared_ptr<bond::ss const> bond)
{ push_sort(_sss, bond); }

void pairbond::push(shared_ptr<bond::vdw const> bond)
{ push_sort(_vdws, bond); }

void pairbond::push(shared_ptr<bond::ionic const> bond)
{ push_sort(_ionics, bond); }

void pairbond::push(shared_ptr<bond::hydrogen const> bond)
{ push_sort(_hydrogens, bond); }

void pairbond::push(shared_ptr<bond::pication const> bond)
{ push_sort(_pications, bond); }

void pairbond::push(shared_ptr<bond::pipistack const> bond)
{ push_sort(_pipistacks, bond); }

void pairbond::push(shared_ptr<bond::generico const> bond)
{ push_sort(_generics, bond); }

pairbond& network::find(chemical_entity::aminoacid const& a, chemical_entity::aminoacid const& b)
{
    auto const key = prelude::concat_lexicographically(a.id(), b.id());
    if (pairbonds_map.find(key) == pairbonds_map.end())
    {
        auto p = std::make_shared<pairbond>();
        pairbonds_map.insert({key, p});
        return *p;
    }

    else
        return *pairbonds_map[key];
}

std::list<shared_ptr<bond::base const>> network::get_one() const
{
    std::list<shared_ptr<bond::base const>> tmp;
    for (auto const& kv: pairbonds_map)
    {
        auto bond = kv.second->get_one();
        if (bond != nullptr)
            tmp.push_back(bond);
    }

    return tmp;
}

std::list<shared_ptr<bond::base const>> network::get_all() const
{
    std::list<shared_ptr<bond::base const>> tmp;
    for (auto const& kv: pairbonds_map)
    {
        auto lst = kv.second->get_all();
        tmp.splice(tmp.end(), lst);
    }

    return tmp;
}

std::list<shared_ptr<bond::base const>> network::get_multiple() const
{
    std::list<shared_ptr<bond::base const>> tmp;
    for (auto const& kv: pairbonds_map)
    {
        auto lst = kv.second->get_multiple();
        tmp.splice(tmp.end(), lst);
    }

    return tmp;
}


std::list<shared_ptr<bond::base const>> pairbond::get_all() const
{
    std::list<shared_ptr<bond::base const>> tmp;

    tmp.insert(tmp.end(), _hydrogens.begin(), _hydrogens.end());
    tmp.insert(tmp.end(), _sss.begin(), _sss.end());
    tmp.insert(tmp.end(), _vdws.begin(), _vdws.end());
    tmp.insert(tmp.end(), _pications.begin(), _pications.end());
    tmp.insert(tmp.end(), _pipistacks.begin(), _pipistacks.end());
    tmp.insert(tmp.end(), _ionics.begin(), _ionics.end());
    tmp.insert(tmp.end(), _generics.begin(), _generics.end());

    return tmp;
}

std::list<shared_ptr<bond::base const>> pairbond::get_multiple() const
{
    std::list<shared_ptr<bond::base const>> tmp;

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
shared_ptr<bond::base const> __best(std::list<shared_ptr<Bond const>> const& l, shared_ptr<bond::base const> b)
{
    static_assert(std::is_base_of<bond::base, Bond>::value, "template typename Bond must inherit from type bond::base");
    if (!l.empty() && (b == nullptr || *b > *l.front()))
    {
        return l.front();
    }
    else
    {
        return b;
    }
}

shared_ptr<bond::base const> pairbond::get_one() const
{
    shared_ptr<bond::base const> best = nullptr;

    best = __best(_hydrogens, best);
    best = __best(_sss, best);
    best = __best(_vdws, best);
    best = __best(_ionics, best);
    best = __best(_pications, best);
    best = __best(_pipistacks, best);
    best = __best(_generics, best);

    return best;
}
