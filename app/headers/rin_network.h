#pragma once

#include <list>
#include <limits>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <type_traits>

namespace chemical_entity
{
class aminoacid;
}

namespace bonds
{
class base;

class hydrogen;

class vdw;

class ss;

class ionic;

class pipistack;

class pication;

class generico;
}

class pairbond final
{
private:
    std::list<bonds::ss const*> _sss;
    std::list<bonds::vdw const*> _vdws;
    std::list<bonds::ionic const*> _ionics;
    std::list<bonds::hydrogen const*> _hydrogens;
    std::list<bonds::pication const*> _pications;
    std::list<bonds::pipistack const*> _pipistacks;
    std::list<bonds::generico const*> _generics;

public:
    ~pairbond();

    void push(bonds::ss const& bond);

    void push(bonds::vdw const& bond);

    void push(bonds::ionic const& bond);

    void push(bonds::hydrogen const& bond);

    void push(bonds::pication const& bond);

    void push(bonds::pipistack const& bond);

    void push(bonds::generico const& bond);

    [[nodiscard]]
    std::list<bonds::base const*> get_multiple() const;

    [[nodiscard]]
    std::list<bonds::hydrogen const*> get_hydrogens() const
    { return _hydrogens; }

    [[nodiscard]]
    std::list<bonds::base const*> get_all() const;

    [[nodiscard]]
    bonds::base const* get_one() const;

    [[nodiscard]]
    bool has_vdw()
    { return !_vdws.empty(); }

    [[nodiscard]]
    bool has_pipi()
    { return !_pipistacks.empty(); }

    [[nodiscard]]
    bool has_backbone()
    { return !_generics.empty(); }
};

class network final
{
private:
    std::unordered_map<std::string, pairbond*> pairbonds_map;

public:
    pairbond& find(chemical_entity::aminoacid const& a, chemical_entity::aminoacid const& b);

    ~network();

    /*
    template<typename Bond, typename... Args>
    void new_bond(Args&& ... args)
    {
        static_assert(std::is_base_of<bonds::base, Bond>::value, "template typename Bond must inherit from type bond");
        auto const* bond = new Bond(args...);

        find(bond->id()).push(*bond);
    }
    */

    std::list<bonds::base const*> get_one() const;

    std::list<bonds::base const*> get_all() const;

    std::list<bonds::base const*> get_multiple() const;

    std::list<bonds::base const*> filter_hbond_realistic(std::list<bonds::base const*> const& input) const;
};