#include <iostream>
#include "error.h"

using namespace boost::python;
using namespace boost;

void handle_pyerror()
{

    PyObject *ptype, *pvalue, *ptraceback;
    PyErr_Fetch(&ptype, &pvalue, &ptraceback);

    handle<> hType(ptype);
    object extype(hType);
    handle<> hTraceback(ptraceback);
    object traceback(hTraceback);

    //Extract error message
    string strErrorMessage = extract<string>(pvalue);
    std::cerr << strErrorMessage<<std::endl;

    //Extract line number (top entry of call stack)
    // if you want to extract another levels of call stack
    // also process traceback.attr("tb_next") recurently
    string filename = extract<string>(traceback.attr("tb_frame").attr("f_code").attr("co_filename"));
    string funcname = extract<string>(traceback.attr("tb_frame").attr("f_code").attr("co_name"));
    std::cerr <<filename<<std::endl;
    std::cerr <<funcname<<std::endl;

return ;
}
