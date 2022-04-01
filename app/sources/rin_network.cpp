#include "rin_network.h"

std::list<bonds::base const*> network::pairbond::get_all() const
{
    std::list<bonds::base const*> tmp;

    tmp.insert(tmp.end(), _hydrogens.begin(), _hydrogens.end());
    tmp.insert(tmp.end(), _sss.begin(), _sss.end());
    tmp.insert(tmp.end(), _vdws.begin(), _vdws.end());
    tmp.insert(tmp.end(), _pications.begin(), _pications.end());
    tmp.insert(tmp.end(), _pipistacks.begin(), _pipistacks.end());
    tmp.insert(tmp.end(), _ionics.begin(), _ionics.end());
    tmp.insert(tmp.end(), _generics.begin(), _generics.end());

    return tmp;
}

std::list<bonds::base const*> network::pairbond::get_multiple() const
{
    std::list<bonds::base const*> tmp;

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
bonds::base const* __best(std::list<const Bond*> const& l, bonds::base const* b)
{
    static_assert(std::is_base_of<bonds::base, Bond>::value, "template typename Bond must inherit from type bond");
    if (!l.empty() && (b == nullptr || *b > *l.front()))
    {
        return l.front();
    }
    else
    {
        return b;
    }
}

bonds::base const* network::pairbond::get_one() const
{
    bonds::base const* best = nullptr;

    best = __best(_hydrogens, best);
    best = __best(_sss, best);
    best = __best(_vdws, best);
    best = __best(_ionics, best);
    best = __best(_pications, best);
    best = __best(_pipistacks, best);
    best = __best(_generics, best);

    return best;
}