import Py4ti2int64_cc
from Py4ti2int64_cc import *

# from Sebastian Gustche PyNormaliz Cone.init() function
def process_function_args(args, kwargs): 
    input_list = [k for k in args]
    for i in kwargs:
        current_input = kwargs[i]
        if type(current_input) == list and len(current_input) > 0 and type(current_input[0]) != list:
            kwargs[i] = current_input
        elif type(current_input) == bool and current_input == True:
            kwargs[i] = current_input = [[]]
        elif type(current_input) == bool and current_input == False:
            kwargs.pop(i)
    for i in kwargs:
        input_list.append(i)
        input_list.append(kwargs[i])

    return input_list

def minimize(*args, **kwargs):
    arguments_list = process_function_args(args, kwargs)
    return Py4ti2int64_cc._minimize(*arguments_list)

def groebner(*args, **kwargs):
    arguments_list = process_function_args(args, kwargs)
    return Py4ti2int64_cc._groebner(*arguments_list)

def normalform(*args, **kwargs):
    arguments_list = process_function_args(args, kwargs)
    return Py4ti2int64_cc._normalform(*arguments_list)

def markov(*args, **kwargs):
    arguments_list = process_function_args(args, kwargs)
    return Py4ti2int64_cc._markov(*arguments_list)

def zbasis(*args, **kwargs):
    arguments_list = process_function_args(args, kwargs)
    return Py4ti2int64_cc._zbasis(*arguments_list)

def walk(*args, **kwargs):
    arguments_list = process_function_args(args, kwargs)
    return Py4ti2int64_cc._walk(*arguments_list)

def hilbert(*args, **kwargs):
    arguments_list = process_function_args(args, kwargs)
    return Py4ti2int64_cc._hilbert(*arguments_list)

def graver(*args, **kwargs):
    arguments_list = process_function_args(args, kwargs)
    return Py4ti2int64_cc._graver(*arguments_list)

def zsolve(*args, **kwargs):
    arguments_list = process_function_args(args, kwargs)
    return Py4ti2int64_cc._zsolve(*arguments_list)

def circuits(*args, **kwargs):
    arguments_list = process_function_args(args, kwargs)
    return Py4ti2int64_cc._circuits(*arguments_list)

def qsolve(*args, **kwargs):
    arguments_list = process_function_args(args, kwargs)
    return Py4ti2int64_cc._qsolve(*arguments_list)

def rays(*args, **kwargs):
    arguments_list = process_function_args(args, kwargs)
    return Py4ti2int64_cc._rays(*arguments_list)
