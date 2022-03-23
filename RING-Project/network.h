#pragma once

#include <list>
#include <limits>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <type_traits>

#include "bonds.h"
#include "parameters.h"

// classe che gestisce l'intera residue interaction network
//
class network
{
private:
	// contiene (informazioni su) tutti i legami tra una coppia di amminoacidi
	// li divide per tipo di legame
	//
	class pairbond
	{
	private:
		// inserisce il legame in fronte alla lista se è migliore del legame in fronte,
		// in fondo altrimenti.
		//
		template <typename Bond>
		static void __push(std::list<const Bond*>& l, Bond const& b)
		{
			static_assert(std::is_base_of<bonds::base, Bond>::value, "template typename Bond must inherit from type bond");
			if (l.empty() || b < *l.front())
				l.push_front(&b);
			else
				l.push_back(&b);
		}

	private:
		std::list<bonds::ss const*> _sss;
		std::list<bonds::vdw const*> _vdws;
		std::list<bonds::ionic const*> _ionics;
		std::list<bonds::hydrogen const*> _hydrogens;
		std::list<bonds::pication const*> _pications;
		std::list<bonds::pipistack const*> _pipistacks;
		std::list<bonds::generico const*> _generics;

	public:
		void push(bonds::ss const& bond) { __push(_sss, bond); }
		void push(bonds::vdw const& bond) { __push(_vdws, bond); }
		void push(bonds::ionic const& bond) { __push(_ionics, bond); }
		void push(bonds::hydrogen const& bond) { __push(_hydrogens, bond); }
		void push(bonds::pication const& bond) { __push(_pications, bond); }
		void push(bonds::pipistack const& bond) { __push(_pipistacks, bond); }
		void push(bonds::generico const& bond) { __push(_generics, bond); }

		// ottiene i legami migliori per tipo
		//
		std::list<bonds::base const*> get_multiple() const;

		std::list<bonds::hydrogen const*> get_hydrogens()
		{
			return _hydrogens;
		}

		// ottiene tutti i legami
		//
		std::list<bonds::base const*> get_all() const;

		// ottiene il miglior legame
		//
		bonds::base const* get_one() const;

	public:
		~pairbond()
		{
			for (bonds::base const* bond : _hydrogens)
				delete bond;

			for (bonds::base const* bond : _sss)
				delete bond;

			for (bonds::base const* bond : _vdws)
				delete bond;

			for (bonds::base const* bond : _pications)
				delete bond;

			for (bonds::base const* bond : _pipistacks)
				delete bond;

			for (bonds::base const* bond : _ionics)
				delete bond;

			for (bonds::base const* bond : _generics)
				delete bond;
		}
	};

private:
	// sostanzialmente è la rin
	//
	std::unordered_map<std::string, pairbond*> pairbonds_map;

	// ritorna il pairbond associato alla chiave corrispondente
	// se non è presente => lo crea, lo inserisce e lo ritorna
	//
	pairbond& select_pairbond(std::string const& get_id)
	{
		if (pairbonds_map.find(get_id) == pairbonds_map.end())
		{
			auto* p = new pairbond();
			pairbonds_map.insert({ get_id, p });
			return *p;
		}

		else
			return *pairbonds_map[get_id];
	}

public:
	~network() { for (auto kv : pairbonds_map) delete kv.second; }

public:
	// crea un nuovo legame e lo inserisce nel pairbond corrispondente alla coppia di amminoacidi che lega
	//
	template <typename Bond, typename... Args>
	void new_bond(Args&&... args)
	{
		static_assert(std::is_base_of<bonds::base, Bond>::value, "template typename Bond must inherit from type bond");
		Bond const* bond = new Bond(args...);
		select_pairbond(bond->id()).push(*bond);
		// verbosity: qui possiamo pretty-printare il bond...
	}

public:
	// ottiene il legame migliore di ogni pairbond della rin
	//
	std::list<bonds::base const*> get_one() const
	{
		std::list<bonds::base const*> tmp;
		for (auto kv : pairbonds_map)
		{
			auto* bond = kv.second->get_one();
			if (bond != nullptr)
				tmp.push_back(bond);
		}

		return tmp;
	}

	// ottiene tutti i legami della rin
	//
	std::list<bonds::base const*> get_all() const
	{
		std::list<bonds::base const*> tmp;
		for (auto kv : pairbonds_map)
		{
			auto lst = kv.second->get_all();
			tmp.splice(tmp.end(), lst);
		}

		return tmp;
	}

	// ottiene i legami migliori per ogni tipo in ogni pairbond della rin
	//
	std::list<bonds::base const*> get_multiple() const
	{
		std::list<bonds::base const*> tmp;
		for (auto kv : pairbonds_map)
		{
			auto lst = kv.second->get_multiple();
			tmp.splice(tmp.end(), lst);
		}

		return tmp;
	}

	std::list<bonds::base const*> filter_hbond_realistic(list<bonds::base const*> const& input) const
	{
		std::set<bonds::hydrogen const*> hydrogen_bonds_output;
		std::vector<bonds::hydrogen const*> hydrogen_bonds_input;
		std::unordered_map<entities::atom const*, int> donors_bond_count;
		std::unordered_map<entities::atom const*, int> hydrogen_bond_count;
		std::unordered_map<entities::atom const*, int> acceptors_bond_count;

		//Get the bonds count of an atom
		auto get_bond_count = [](std::unordered_map<entities::atom const*, int> &container, entities::atom const* atom) -> int
		{
			if (container.find(atom) == container.end())
				return 0;
			else
				return container[atom];
		};
		//Increase the bonds count of an atom
		auto inc_bond_count = [](std::unordered_map<entities::atom const*, int> &container, entities::atom const* atom) -> void
		{
			if (container.find(atom) == container.end())
				container[atom] = 0;
			container[atom]++;
		};
		auto can_be_added = [&](bonds::hydrogen const* bond) -> bool
		{
			return (get_bond_count(donors_bond_count, bond->donor_ptr()) < bond->donor().how_many_hydrogen_can_donate() &&
					get_bond_count(hydrogen_bond_count, bond->hydrogen_ptr()) < 1 && //An hydrogen can make only one bond
					get_bond_count(acceptors_bond_count, bond->acceptor_ptr()) < bond->acceptor().how_many_hydrogen_can_accept());
		};
		auto add_bond = [&](bonds::hydrogen const* bond) -> void
		{
			inc_bond_count(donors_bond_count, bond->donor_ptr());
			inc_bond_count(hydrogen_bond_count, bond->hydrogen_ptr());
			inc_bond_count(acceptors_bond_count, bond->acceptor_ptr());
			hydrogen_bonds_output.insert(bond);
		};

		//Extract hydrogen bonds from input
		for (auto& i : input)
		{
			if(i->get_type() == "hydrogen")
				hydrogen_bonds_input.push_back((bonds::hydrogen const*)i);
		}

		//Order from smallest to largest energy
		sort(hydrogen_bonds_input.begin(), hydrogen_bonds_input.end(), [](bonds::hydrogen const* a, bonds::hydrogen const* b) { return a->get_energy() < b->get_energy(); });
		
		//Add as many hydrogen bonds as possible
		for (auto i : hydrogen_bonds_input)
		{
			if (can_be_added(i))
				add_bond(i);
		}

		//Let's build the output list
		std::list<bonds::base const*> output;
		for (auto i : input)
		{
			//Insert i into the output if it is not an hydrogen or if it is in the filtered list
			if (i->get_type() != "hydrogen" || hydrogen_bonds_output.find((bonds::hydrogen const*)i) != hydrogen_bonds_output.end())
			{
				output.push_back(i);
			}
		}
		return output;
	}
};