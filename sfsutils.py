"""Parser for KSP save file"""
from collections import OrderedDict
import re

def parse_savefile(sfs, sfs_is_path=True):
    """Parses an SFS file
    Params:
        sfs: str; the path to the SFS file to read or a string containing data read from an sfs.
        sfs_is_path (optional, default True): bool; whether the 'sfs' param is a path or raw data.
    Raises:
        No specific exceptions.
    Returns:
        OrderedDict containing the data in the SFS.
    Extra information:
        All values are strings as SFS files do not reveal data to be any type.
        The SFS format is particularly bad and this leads to the returned OrderedDict
        containing data that is unusually structured. If the SFS contains multiple keys of any
        kind with the same name (this can be a 'node' header or values in a node), then the data
        contained within these keys will formatted as the common name of the keys as a key
        in a dict, and the values as a list. This data will always be in the exact order
        that they were in in the SFS. Example:
        --SFS format--
        NODE
        {
            x = 1
            x = 2
            y = 3
        }
        NODE
        {
            value = 1
        }
        OTHER
        {
            z = 4
        }
        --Python structure--
        {
            "NODE": [
                {"x": ["1","2"], "y": "3"},
                {"value": "1"}
            ],
            "OTHER": {
                "z": "4"
            }
        }
    """
    if sfs_is_path:
        data = open(sfs, "r").read()
    else:
        data = sfs
    # removes double slash "//" comments
    data = re.sub('(?m)//.*', '', data)
    # removes all tabs and cursor marks
    data = data.replace("    ", "\t")   # remove large numbers of spaces
    data = data.replace("=", " = ")     # ensure there is a space before each =
    data = data.replace("\t", "")       # remove tabs
    data = data.replace("\r", "")       # remove cursor mark
    # removes % chars (for Eeloo in OPM)
    data = data.replace("%", "")
    # replace multiple newlines or spaces with singles
    while "\n\n" in data:
        data = data.replace("\n\n", "\n")
    while "  " in data:
        data = data.replace("  ", " ")
    # in_nodes tracks the location of data being parsed (what nodes the parser is inside)
    in_nodes = []
    out_dict = OrderedDict()
    # key_read contains the start and end index of the key being read
    key_read = [0, None]
    # value_read contains the start index of the value being read
    value_read = None
    trigger = set(("\n", "}", "="))
    key_ignore = set(("\n"," ", "{", "}", "="))
    for index, char in enumerate(data[0:-1]):
        # check if the char is one of the chars which leads to an action
        # this is an optimisation only
        if char in trigger:
            if char == "\n":
                # if the key is empty, continue
                if (key_read[0] == index - 1) and (data[key_read[0]] in key_ignore):
                    pass
                # if next char is an open bracket, save it as a new node
                else:
                    if data[index + 1] == "{":
                        in_nodes.append(data[key_read[0]: index])
                        write_list = in_nodes[:]
                        write_list.append(OrderedDict())
                    # else it is a value in an existing node
                    else:
                        write_list = in_nodes[:]
                        write_list.append(data[key_read[0]: key_read[1]])
                        write_list.append(data[value_read: index])
                    set_value(out_dict, write_list)
                if data[index-1] in key_ignore:
                    key_read = [index + 1, None]
                else:
                    key_read = [index + 1, None]
            # pop the end of the 'stack' used to track attribute location
            # when the end of a node is found
            elif char == "}":
                try: 
                    if 'vessels2' in in_nodes and 'computer' in in_nodes:
                        print("")
                    in_nodes.pop()
                except IndexError:
                    pass
            # set the end index of the key and start index of the value
            # due to the set check (with 'trigger'), the char must be =
            else:
                key_read[1] = index - 1
                value_read = index + 2
    return out_dict

def set_value(dict_nested, address_list):
    """Sets a value in a nested dict
    WARNING - modifies the dictionary passed as an arg"""
    # references the main dict
    current = dict_nested
    # locate the desired node to write to through iterating through the keys
    # while selecting the last element of any list found, as the data is in order
    for path_item in address_list[:-2]:
        if isinstance(current, list):
            current = current[-1][path_item]
        else:
            current = current[path_item]
    # if current is a list, then take the last entry as that's what will be modified
    if isinstance(current, list):
        current = current[-1]
    # if the node already exists
    if address_list[-2] in current:
        # if it's a list simply append it to the list
        if isinstance(current[address_list[-2]], list):
            current[address_list[-2]].append(address_list[-1])
        # else convert the existing dict to a list
        else:
            current[address_list[-2]] = [current[address_list[-2]], address_list[-1]]
    # if it doesn't exist
    else:
        # guaranteed to be a dict thanks to earlier list check, so insert the key into the dict
        current[address_list[-2]] = address_list[-1]


def writeout_savefile(parsed_data, destination_file=""):
    """Writes out the parsed data back into the SFS format
    Params:
        parsed_data: str; the parsed dictionary generated by parse_savefile.
        destination_file (optional): str; the destination file to write the SFS to.
    Raises:
        No specific exceptions
    Returns:
        str containg the generated SFS if a destination file is not specified
        None if a destination file is specified
    Extra information:
        This function will generate a byte perfect copy of the original SFS parsed assuming
        the data is not modified. All abnormalities of the SFS format are addressed and
        represented correctly.
    """
    indents = -1
    out_str = serialise_data(parsed_data, indents)
    if not destination_file:
        return out_str
    open(destination_file, "w").write(out_str)
    return None

def serialise_data(obj, indents, outer_key=None):
    """Recursively serialises data"""
    # indent upon each recurse
    indents += 1
    out_str = ""
    # set up the buffer list
    buffer_list = []
    if isinstance(obj, list):
        for item in obj:
            buffer_list.append(item)
    else:
        buffer_list.append(obj)
    for item in buffer_list:
        # if it is a string, it is one of SFS stupid same keys
        # to different values, so just write value to node
        if isinstance(item, str):
            out_str += write_value_to_node(indents, outer_key, item)
        else:
            # it is a dict, so iterate through
            for key, value in item.items():
                # if value is a string, it must be a value to write to a node
                if isinstance(value, str):
                    out_str += write_value_to_node(indents, key, value)
                # if it's a dict, it's another node, so recurse
                elif isinstance(value, dict):
                    out_str += write_new_node(indents, key, value)
                # if it's a list it could be multiple things
                elif isinstance(value, list):
                    # if everything in the list is a string, then it is one of the multi
                    # value nodes (could this be optimised TODO)
                    if all(isinstance(x, str) for x in value):
                        out_str += serialise_data(value, indents - 1, outer_key=key)
                    # else just write out each subdict in the list
                    else:
                        for subdict in value:
                            out_str += write_new_node(indents, key, subdict)
    return out_str

def write_new_node(indents, sect_name, value):
    """Write a new node to the SFS"""
    # adds the header
    out_str = "{0}{1}\n{0}{{\n".format("\t" * indents, sect_name)
    # adds data through recursion
    out_str += serialise_data(value, indents, outer_key=sect_name)
    # closes the block
    out_str += "{0}}}\n".format("\t" * indents)
    return out_str

def write_value_to_node(indents, key, value):
    """Writes a key value pair into a node"""
    return "{0}{1} = {2}\n".format("\t" * indents, key, value)