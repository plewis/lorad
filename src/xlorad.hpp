#pragma once

#include <boost/format.hpp>

namespace lorad {

    class XLorad : public std::exception {
        public:
                                XLorad() throw() {}
                                XLorad(const std::string s) throw() : _msg() {_msg = s;}
                                XLorad(const boost::format & f) throw() : _msg() {_msg = boost::str(f);}
            virtual             ~XLorad() throw() {}
            const char *        what() const throw() {return _msg.c_str();}

        private:

            std::string         _msg;
    };

}
