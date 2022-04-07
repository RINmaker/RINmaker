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

namespace bond
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
    std::list<bond::ss const*> _sss;
    std::list<bond::vdw const*> _vdws;
    std::list<bond::ionic const*> _ionics;
    std::list<bond::hydrogen const*> _hydrogens;
    std::list<bond::pication const*> _pications;
    std::list<bond::pipistack const*> _pipistacks;
    std::list<bond::generico const*> _generics;

public:
    ~pairbond();

    void push(bond::ss const& bond);

    void push(bond::vdw const& bond);

    void push(bond::ionic const& bond);

    void push(bond::hydrogen const& bond);

    void push(bond::pication const& bond);

    void push(bond::pipistack const& bond);

    void push(bond::generico const& bond);

    [[nodiscard]]
    std::list<bond::base const*> get_multiple() const;

    [[nodiscard]]
    std::list<bond::hydrogen const*> get_hydrogens() const
    { return _hydrogens; }

    [[nodiscard]]
    std::list<bond::base const*> get_all() const;

    [[nodiscard]]
    bond::base const* get_one() const;

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

    std::list<bond::base const*> get_one() const;

    std::list<bond::base const*> get_all() const;

    std::list<bond::base const*> get_multiple() const;

    std::list<bond::base const*> filter_hbond_realistic(std::list<bond::base const*> const& input) const;
};