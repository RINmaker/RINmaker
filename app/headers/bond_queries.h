#pragma once

#include <array>
#include <string>
#include <set>
#include <tuple>

#include "log_manager.h"

namespace chemical_entity
{
    class aminoacid;
    class atom;
    class ring;
    class ionic_group;
}

class network;

namespace bondfunctors {

class base {
protected:
    network& _net;
    int _nbonds;

protected:
    explicit base(network& net)
            : _net(net), _nbonds(0) {}

public:
    virtual ~base() = default;
};

/**
 * Van der Waals bond test.
 */
class vdw : public base {
private:
    // just to keep track of already accepted pairs; vdw is tested for ALL VDWS against ALL VDWS
    std::set<std::string> _bonded;

public:
    explicit vdw(network& net)
            : base(net) {}

    ~vdw() override { log_manager::main()->info("found {} vdw bonds", _nbonds); }

public:
    void operator()(chemical_entity::atom const& a, chemical_entity::atom const& b);
};

/**
 * Ionic bond test.
 */
class ionic : public base {
public:
    explicit ionic(network& net)
            : base(net) {}

    ~ionic() override { log_manager::main()->info("found {} ionic bonds", _nbonds); }

public:
    // aminoacids must be separated and have opposite charges
    void operator()(chemical_entity::ionic_group const& a, chemical_entity::ionic_group const& b);
};

/**
 * Hydrogen bond test
 */
class hydrogen : public base {
public:
    explicit hydrogen(network& net)
            : base(net) {}

    ~hydrogen() override { log_manager::main()->info("found {} hydrogen bonds", _nbonds); }

public:
    // a bit complex, please refer to the paper
    void operator()(chemical_entity::atom const& acceptor, chemical_entity::atom const& donor);
};

/**
 * PiCation bond test.
 */
class pication : public base {
public:
    explicit pication(network& net)
            : base(net) {}

    ~pication() override { log_manager::main()->info("found {} pication bonds", _nbonds); }

public:
    void operator()(chemical_entity::atom const& cation, chemical_entity::ring const& ring);
};

/**
 * PiPiStacking bond test.
 */
class pipistack : public base {
private:
    // just to keep track of already accepted pairs; pipistack is tested for ALL RINGS against ALL RINGS
    std::set<std::string> _bonded;

public:
    explicit pipistack(network& net)
            : base(net) {}

    ~pipistack() override { log_manager::main()->info("found {} pipistack bonds", _nbonds); }

public:
    void operator()(chemical_entity::ring const& a, chemical_entity::ring const& b);
};

/**
 * Generic bond test (alpha-alpha or beta-beta).
 */
class generico : public base {
private:
    // just to keep track of already accepted pairs; generic is tested for ALL CARBONS against ALL CARBONS
    std::set<std::string> _bonded;

public:
    explicit generico(network& net)
            : base(net) {}

    ~generico() override { log_manager::main()->info("found {} generic bonds", _nbonds); }

public:
    // tecnicamente il test è sempre valido, ricordiamoci che questo è il test applicato ai vicini della range search!
    //
    void operator()(chemical_entity::atom const& a, chemical_entity::atom const& b);
};
}
