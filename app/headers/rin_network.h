#pragma once

#include <list>
#include <limits>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <type_traits>
#include <memory>

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
    std::list<std::shared_ptr<bond::ss const>> _sss;
    std::list<std::shared_ptr<bond::vdw const>> _vdws;
    std::list<std::shared_ptr<bond::ionic const>> _ionics;
    std::list<std::shared_ptr<bond::hydrogen const>> _hydrogens;
    std::list<std::shared_ptr<bond::pication const>> _pications;
    std::list<std::shared_ptr<bond::pipistack const>> _pipistacks;
    std::list<std::shared_ptr<bond::generico const>> _generics;

public:
    void push(std::shared_ptr<bond::ss const> bond);

    void push(std::shared_ptr<bond::vdw const> bond);

    void push(std::shared_ptr<bond::ionic const> bond);

    void push(std::shared_ptr<bond::hydrogen const> bond);

    void push(std::shared_ptr<bond::pication const> bond);

    void push(std::shared_ptr<bond::pipistack const> bond);

    void push(std::shared_ptr<bond::generico const> bond);

    [[nodiscard]]
    std::list<std::shared_ptr<bond::base const>> get_multiple() const;

    [[nodiscard]]
    std::list<std::shared_ptr<bond::hydrogen const>> get_hydrogens() const
    { return _hydrogens; }

    [[nodiscard]]
    std::list<std::shared_ptr<bond::base const>> get_all() const;

    [[nodiscard]]
    std::shared_ptr<bond::base const> get_one() const;

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
    std::unordered_map<std::string, std::shared_ptr<pairbond>> pairbonds_map;

public:
    pairbond& find(chemical_entity::aminoacid const& a, chemical_entity::aminoacid const& b);

    std::list<std::shared_ptr<bond::base const>> get_one() const;

    std::list<std::shared_ptr<bond::base const>> get_all() const;

    std::list<std::shared_ptr<bond::base const>> get_multiple() const;

};