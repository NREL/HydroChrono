import os



from docutils import nodes
from docutils.nodes import Node, make_id

from docutils.parsers.rst import directives, Directive
from docutils.parsers.rst.directives import images, tables
from docutils.parsers.rst.directives.tables import Table

from sphinx.util import logging

LOG = logging.getLogger(__name__)

import h5py

import pdb


class HydroChrono_H5(Table):

    option_spec = Table.option_spec
    option_spec['h5-file'] =  directives.unchanged_required
    option_spec['name'] =  directives.unchanged

    required_arguments = 1

    def run(self):
        
        env = self.state.document.settings.env
        app = env.app
        config = app.config


        LOG.info('Running ..... %s', self.options.get('h5-file') )

        MAX_COL = 3
        col_widths = self.get_column_widths(MAX_COL)
        title, messages = self.make_title()

        LOG.info('TITLE ..... %s', title )
        LOG.info('MESHG ..... %s', messages )

        table_headers = ['AAA', 'BBB', 'CCC']
        table_data = [[1,2,3], ['quatre', 'cinq', 'six']]


        h5f = h5py.File(self.options.get('h5-file'), 'r')

        dset = h5f['simulation_parameters/T']
        for k in h5f.keys():
            LOG.info('%s', k)

        for k in dset:
            LOG.info('%s', k)


        table_node = self._build_table(table_data, col_widths, table_headers)


        # pdb.set_trace()
        # 
        
        self.add_name(table_node)
        if title: 
            table_node.insert(0, title) 

        return [table_node] + messages
        #return super().run()


    def _build_table(self, table_data, col_widths, headers):
        
        table = nodes.table()

        tgroup = nodes.tgroup(cols=len(col_widths))
        table += tgroup
        tgroup.extend(nodes.colspec(colwidth=col_width) for col_width in col_widths)
        
        # Headers
        thead = nodes.thead()
        tgroup += thead
        row_node = nodes.row()
        thead += row_node
        row_node.extend(
            nodes.entry(h, nodes.paragraph(text=h)) for h in headers
        )

        # Rows
        tbody = nodes.tbody()
        tgroup += tbody
        rows = []
        for row in table_data:
            trow = nodes.row()
            for val in row:
                entry = nodes.entry()
                entry += nodes.paragraph(text=str(val))
                trow += entry
            
            rows.append(trow)

        tbody.extend(rows)

        return table

class HydroChrono_H5Bis(Directive):
    
    has_content = True
    required_arguments = 0
    optional_arguments = 1
    final_argument_whitespace = False

    def run(self):

        document = self.state_machine.document


        setup.app.builder.srcdir
        #directives.uri(arguments[0]))

        srcfile = document.attributes['source']
        rst_file = document.attributes['source']
        rst_dir = os.path.dirname(rst_file)

        text=f'source: {srcfile}   \n \n -- path: {rst_dir}   path2: {setup.app.builder.srcdir}'

        paragraph_node = nodes.paragraph(text=text) # f'bla bla lllllll {rst_file} \n {self.arguments[0]}')
        paragraph_node2 = nodes.paragraph(text=text)

        return [paragraph_node, paragraph_node2]

        #try:
            #rst_file = self.state_machine.document.attributes['source']
            #rst_dir = os.path.dirname(rst_file)

            #self.arguments
            #self.content
            #self.options
            #self.state_machine
            #self.state
            #self.lineno

            #paragraph_node = nodes.paragraph(text='bibibibib') # f'bla bla lllllll {rst_file} \n {self.arguments[0]}')

            #return [paragraph_node]
        
        #except Exception as e:
        #    raise self.error(str(e))



def setup(app):

    setup.app = app

    app.add_directive("hydrochrono_h5", HydroChrono_H5)

    LOG.info('Setting .....')
    
    return {

        'version': '0.1',

        'parallel_read_safe': True,

        'parallel_write_safe': True,

    }


list1 = [
    "simulation_parameters/rho",
    "simulation_parameters/g",
    "simulation_parameters/water_depth",
    "simulation_parameters/w"
]

liste2 = [
    "/properties/disp_vol",
    "/hydro_coeffs/radiation_damping/impulse_response_fun/t",
    "/properties/cg",
    "/properties/cb",
    "/hydro_coeffs/linear_restoring_stiffness",
    "/hydro_coeffs/added_mass/inf_freq",
    "/hydro_coeffs/radiation_damping/impulse_response_fun/K",
    "/hydro_coeffs/excitation/mag",
    "/hydro_coeffs/excitation/phase",
    "/hydro_coeffs/excitation/impulse_response_fun/t",
    "/hydro_coeffs/excitation/impulse_response_fun/f"
]

def _h5_body_list(h5f):
    import re
    rebody = re.compile("body(\d+)")

    for key in h5f.keys():
        m = rebody.match(key)
        if m:    
            print(f"found '{key}'  with group '{m.group(1)}' mame: '{h5f[key+'/properties/name'].asstr()[...]}'")

    return []

def _h5hydo_read(filename: str):

    h5f = h5py.File(filename, 'r')

    # BEM Data (Wamit or other types ...)
    dset_bem = h5f['bem_data/code']
    print(dset_bem.name)
    print(dset_bem.asstr()[...])
    print(dset_bem.shape)
    for att in dset_bem.attrs:
        print(att)    


    dset = h5f['simulation_parameters/T']
    for k in h5f.keys():
        print(f'{k}')

    #for k in dset:
    #    print(f'{k}')


    _h5_body_list(h5f)


if __name__ == '__main__':
    import sys
    print('Start in main')
    _h5hydo_read(sys.argv[1])

