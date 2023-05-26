import os

import numpy as np
import pandas as pd

from lxml import objectify
from lxml.etree import tostring, Element, SubElement


class DataHandler:

    @staticmethod
    def read_file(filename):
        '''Opens a text file, returns the lines as list.'''
        with open(filename, 'r') as f:
            lines = f.read().splitlines()
        return lines

    @staticmethod
    def trim_XML_string(lines):
        '''Remove declaration line and remove new-line and indentation signs. Joins the string list into one string.'''
        lines = lines[1:]
        splitlines = [line.replace('\t', '') for line in lines]
        XML_trimmed = ''.join(splitlines)
        return XML_trimmed

    @staticmethod
    def create_objectify_from_XML_file(filename):
        '''Reads the .xml file and creates an lxml.objectify.ObjectifiedElement.'''
        lines = DataHandler.read_file(filename)
        XML_trimmed = DataHandler.trim_XML_string(lines)
        obj = objectify.fromstring(XML_trimmed)
        return obj

    @staticmethod
    def getXpath(element):
        path = ''
        ancestor_lst = []
        for ancestor in element.iterancestors():
            ind = ''
            if ancestor_lst:
                if len(ancestor.findall(ancestor_lst[-1].tag)) > 1:
                    ind = f'[{ancestor.findall(ancestor_lst[-1].tag).index(ancestor_lst[-1])}]'
            attrib = ' and '.join([f"@{k}='{v}'" for k, v in ancestor.attrib.items()])
            if attrib:
                path = '/' + ancestor.tag + f"[{attrib}]" + ind + path
            else:
                path = '/' + ancestor.tag + ind + path
            ancestor_lst.append(ancestor)
        return path + '/' + element.tag

    @staticmethod
    def get_variables_in_dataschema_object(obj):
        def getVariables(obj, oldpath='', depth=0, variables=[]):
            '''Recursive function that returns all variables (lowest level elements that contain values) in a list.
            Returns of every element:
                -tag
                -xpath
                -base value
            '''
            depth += 1
            for i in obj.iterchildren():
                if type(i) != objectify.ObjectifiedElement:
                    path = DataHandler.getXpath(i).replace(i.tag, "")
                    if path != oldpath:
                        oldpath = path
                    variables.append([i.tag, oldpath + i.tag, i.text])

                variables = getVariables(i, oldpath, depth, variables)
            return variables

        variables_in_ds = getVariables(obj)

        return variables_in_ds

    @staticmethod
    def add_subelement(parent_ele, tag, value=None, XPath=None):
        child_ele = SubElement(parent_ele, tag)
        if value:
            child_ele.text = str(value)
        if XPath:
            child_ele.set('XPath', XPath)
        return child_ele

    @staticmethod
    def add_variable_element(parent_ele, tag, values, XPath):
        child_ele = DataHandler.add_subelement(parent_ele=parent_ele, tag=tag, XPath=XPath)
        DataHandler.add_subelement(parent_ele=child_ele, tag='Values', value=values)
        DataHandler.add_subelement(parent_ele=child_ele, tag='Lower_bound', value=min(values))
        DataHandler.add_subelement(parent_ele=child_ele, tag='Upper_bound', value=max(values))
        return child_ele

    @staticmethod
    def combine_instances_to_dataset(foldername, outputfilename, inputs=None):
        '''
        Function to combine the multiple output XMLs produced by RCE into one file.

        In case a custom design table is chosen, RCE does not seem to map the input values to the output XMLs.
        Thus, they have to be inserted retroactively.

        Since this function is designed to be used in combination with the RCEInterface of the SAS module,
        the inputs follow a similar format and can thus easily be reused.

        :param foldername: Name of or path to the folder that contains the results produced by RCE.
        :param outputfilename: Filename of the output, the combined XML file.
        :param inputs: Inputs of the DOE and their samples. Format {'Variable name': [sample1, sample2, ..., sampleN]}
        :return:
        '''

        print(f'Retrieving storage folder {foldername}')

        schema_objs = []

        for i, filename in enumerate(os.listdir(foldername)):
            obj = DataHandler.create_objectify_from_XML_file(f'{foldername}\{filename}')
            schema_objs.append(obj)

        print(f'Found {(n_schema_objs := len(schema_objs))} files to combine')

        # for obj in schema_objs:
        #     print(obj.aircraft.geometry.find('lambda').text)

        '''Determine if parameters are consistent across the entire dataset'''
        varsdict = []

        for i, obj in enumerate(schema_objs):
            variables_of_ds = DataHandler.get_variables_in_dataschema_object(obj)
            vardict = {}
            for var in variables_of_ds:
                # print(var)
                vardict.update({
                    var[1]: var[2]
                })
            varsdict.append(vardict)

            '''Create set of XPaths of the first dataschema instance. XPaths sets of the other instances should exactly equal to this. '''
            if i == 0:
                XPaths_in_ds = list(vardict.keys())
            else:
                assert (
                set(vardict.keys()) == set(XPaths_in_ds), 'Data instances should follow the same schema! Aborting...')

        '''Create pandas dataframe with elements and their values'''

        data_df = pd.DataFrame(columns=XPaths_in_ds,
                               index=range(1, n_schema_objs + 1))

        for i, vardict in enumerate(varsdict):
            for xpath, value in vardict.items():
                data_df[xpath][i + 1] = value

        '''Dataframe has been created. Can also be returned for reference?'''
        data_dicts = data_df.to_dict(orient='list')

        '''Determine constants and variables across set'''
        const_dict = {}
        inp_dict = {}
        out_dict = {}

        '''Determine consistency across values'''
        for parameter, values_lst in data_dicts.items():
            if len(set(values_lst)) != 1:
                # print(f'Parameter {line[0]} varies across the dataset with values {val_lst}')
                values_lst_float = [float(i) for i in values_lst]
                if inputs:
                    if parameter in inputs.keys():
                        inp_dict.update({parameter: values_lst_float})
                    else:
                        out_dict.update({parameter: values_lst_float})
                else:
                    inp_dict.update({parameter: values_lst_float})

            else:
                # print(f'Parameter {line[0]} is constant across the dataset with value {set(val_lst)}')
                const_dict.update({parameter: list(set(values_lst))[0]})

        '''Create data schema object'''
        DSESchema = Element('DOESchema')

        HEADER = 'Header'
        CONSTANT = 'Constant'

        header_dict = dict(
            author='Author here',
            num_designpoints=str(len(schema_objs))
        )

        '''If I/O is known, make distinction'''
        if inputs:
            INPUT = 'Input'
            OUTPUT = 'Output'

            children = [HEADER, CONSTANT, INPUT, OUTPUT]

            '''Repetive code, can be cleaned up with a function '''
            for child in children:
                child_ele = SubElement(DSESchema, child)
                if child == HEADER:
                    for headerkey, value in header_dict.items():
                        DataHandler.add_subelement(child_ele, headerkey, value)
                elif child == CONSTANT:
                    for XPath, value in const_dict.items():
                        const_name = XPath.split('/')[-1]
                        DataHandler.add_subelement(child_ele, const_name, value, XPath)
                elif child == INPUT:
                    for XPath, values in inp_dict.items():
                        inp_name = XPath.split('/')[-1]
                        DataHandler.add_variable_element(child_ele, inp_name, values, XPath)

                elif child == OUTPUT:
                    for XPath, values in out_dict.items():
                        out_name = XPath.split('/')[-1]
                        DataHandler.add_variable_element(child_ele, out_name, values, XPath)

        else:
            VARIABLE = 'Variable'

            children = [HEADER, CONSTANT, VARIABLE]

            for child in children:
                child_ele = SubElement(DSESchema, child)
                if child == HEADER:
                    for headerkey, value in header_dict.items():
                        DataHandler.add_subelement(child_ele, headerkey, value)
                elif child == CONSTANT:
                    for XPath, value in const_dict.items():
                        const_name = XPath.split('/')[-1]
                        DataHandler.add_subelement(child_ele, const_name, value, XPath)
                else:
                    for XPath, values in inp_dict.items():
                        inp_name = XPath.split('/')[-1]
                        DataHandler.add_variable_element(child_ele, inp_name, values, XPath)

        ''''Write to document'''
        filename = f'{outputfilename}.xml'
        with open(filename, 'wb') as doc:
            doc.write(tostring(DSESchema, xml_declaration=True, encoding='utf-8', pretty_print=True))
        return

    @staticmethod
    def search_dataschema_for_basevalues(obj,varnames):
        basevalues = dict()
        for var in varnames:
            found = False
            for cat in obj.iterchildren():
                if cat.find(var) != None:
                    basevalues.update({var: cat.find(var).text})
                    break
                if not found:
                    for child in cat.iterchildren():
                        if child.find(var) != None:
                            basevalues.update({var: child.find(var).text})
                            break
        return basevalues