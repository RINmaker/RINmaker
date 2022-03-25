#pragma once

#include <string>
#include <array>

#include "prelude.h"
#include "utils/interval.h"

namespace records {
    template<size_t N, typename Derived>
    class base {
    protected:
        std::array<std::string, N> _fields;

        explicit base(std::string &line) {
            size_t line_size = line.size();
            for (size_t i = 0; i < N; ++i) {
                size_t first = Derived::_limits[i].first;
                size_t second = Derived::_limits[i].second;

                // TODO: throw
                if (first + second <= line_size)
                    _fields[i] = prelude::trim(line.substr(first, second));
            }
        }
    };

    class atom : public base<15, atom> {
    private:
        friend class base<15, atom>;

        static std::array<std::pair<size_t, size_t>, 15> const _limits;

    public:
        explicit atom(std::string &line) : base(line) {}

        [[nodiscard]]
        bool same_res(atom const &other) const {
            return chain_id() == other.chain_id() && res_name() == other.res_name() && res_seq() == other.res_seq();
        }

        [[nodiscard]]
        std::string const &name() const { return _fields[2]; }

        [[nodiscard]]
        std::string const &res_name() const { return _fields[4]; }

        [[nodiscard]]
        std::string const &chain_id() const { return _fields[5]; }

        [[nodiscard]]
        std::string const &element_name() const { return _fields[13]; }

        [[nodiscard]]
        int serial() const { return std::stoi(_fields[1]); }

        [[nodiscard]]
        int charge() const {
            auto const &x = _fields[14];
            return x == "1+" ? 1 : x == "1-" ? -1 : 0;
        } // @TODO config

        [[nodiscard]]
        int res_seq() const { return std::stoi(_fields[6]); }

        [[nodiscard]]
        double x() const { return std::stod(_fields[8]); }

        [[nodiscard]]
        double y() const { return std::stod(_fields[9]); }

        [[nodiscard]]
        double z() const { return std::stod(_fields[10]); }

        [[nodiscard]]
        double temp_factor() const { return std::stod(_fields[12]); }
    };

    class helix : public base<14, helix> {
    private:
        friend class base<14, helix>;

        static std::array<std::pair<size_t, size_t>, 14> const _limits;

    public:
        explicit helix(std::string &line) : base(line) {}

        [[nodiscard]]
        std::string const &init_chain_id() const { return _fields[4]; }

        [[nodiscard]]
        std::string const &init_res_name() const { return _fields[3]; }

        [[nodiscard]]
        int init_seq_number() const { return std::stoi(_fields[5]); }

        [[nodiscard]]
        std::string const &end_chain_id() const { return _fields[8]; }

        [[nodiscard]]
        std::string const &end_res_name() const { return _fields[7]; }

        [[nodiscard]]
        int end_seq_number() const { return std::stoi(_fields[9]); }

        [[nodiscard]]
        std::string get_id() const { return _fields[2]; }

        [[nodiscard]]
        std::string init_res_id() const {
            return init_chain_id() + ":" + std::to_string(init_seq_number()) + ":_:" + init_res_name();
        }

        [[nodiscard]]
        std::string end_res_id() const {
            return end_chain_id() + ":" + std::to_string(end_seq_number()) + ":_:" + end_res_name();
        }

        [[nodiscard]]
        interval<int> range() const { return {init_seq_number(), end_seq_number()}; }

        [[nodiscard]]
        int serial() const { return std::stoi(_fields[1]); }
    };

    class sheet_piece : public base<23, sheet_piece> {
    private:
        friend class base<23, sheet_piece>;

        static std::array<std::pair<size_t, size_t>, 23> const _limits;

    public:
        explicit sheet_piece(std::string &line) : base(line) {}

        [[nodiscard]]
        std::string const &init_chain_id() const { return _fields[5]; }

        [[nodiscard]]
        std::string const &init_res_name() const { return _fields[4]; }

        [[nodiscard]]
        int init_seq_number() const { return std::stoi(_fields[6]); }

        [[nodiscard]]
        std::string const &end_chain_id() const { return _fields[9]; }

        [[nodiscard]]
        std::string const &end_res_name() const { return _fields[8]; }

        [[nodiscard]]
        int end_seq_number() const { return std::stoi(_fields[10]); }

        [[nodiscard]]
        std::string get_id() const { return _fields[2]; }

        [[nodiscard]]
        std::string init_res_id() const {
            return init_chain_id() + ":" + std::to_string(init_seq_number()) + ":_:" + init_res_name();
        }

        [[nodiscard]]
        std::string end_res_id() const {
            return end_chain_id() + ":" + std::to_string(end_seq_number()) + ":_:" + end_res_name();
        }

        [[nodiscard]]
        interval<int> range() const { return {init_seq_number(), end_seq_number()}; }

        [[nodiscard]]
        int incremental_strand_number() const { return std::stoi(_fields[1]); }
    };

    class ss : public base<13, ss> {
    private:
        friend class base<13, ss>;

        static std::array<std::pair<size_t, size_t>, 13> const _limits;

    public:
        explicit ss(std::string &line) : base(line) {}

        [[nodiscard]]
        double length() const { return std::stod(_fields[12]); }

        [[nodiscard]]
        int seq_num_1() const { return std::stoi(_fields[4]); }

        [[nodiscard]]
        int seq_num_2() const { return std::stoi(_fields[8]); }

        [[nodiscard]]
        std::string const &name_1() const { return _fields[2]; }

        [[nodiscard]]
        std::string const &name_2() const { return _fields[6]; }

        [[nodiscard]]
        std::string const &chain_id_1() const { return _fields[3]; }

        [[nodiscard]]
        std::string const &chain_id_2() const { return _fields[7]; }
    };
}
