# -*- coding: utf-8 -*-


from IPython.core.display import display, Javascript


def init_datatable_mode(pd):
    """Initialize DataTable mode for pandas DataFrame represenation."""

    # configure path to the datatables library using requireJS
    # that way the library will become globally available
    display(Javascript("""
        require.config({
            paths: {
                DT: '//cdn.datatables.net/1.10.19/js/jquery.dataTables.min',
            }
        });
        $('head').append('<link rel="stylesheet" type="text/css" href="//cdn.datatables.net/1.10.19/css/jquery.dataTables.min.css">');
    """))

    def _repr_datatable_(self):
        """Return DataTable representation of pandas DataFrame."""
        # classes for dataframe table (optional)
        classes = ['table', 'table-striped', 'table-bordered']

        # create table DOM
        script = (
            f'$(element).html(`{self.to_html(index=False, classes=classes)}`);\n'
        )

        # execute jQuery to turn table into DataTable
        script += """
            require(["DT"], function(DT) {
                $(document).ready( () => {
                    // Turn existing table into datatable
                    $(element).find("table.dataframe").DataTable();
                })
            });
        """

        return script

    pd.DataFrame._repr_javascript_ = _repr_datatable_


import ipyvuetify as v
from traitlets import (Unicode, List)

def process_index(dataframe_, index_name="gene"):
    dataframe = dataframe_.copy()
    dataframe.index.name = index_name
    dataframe = dataframe.reset_index(drop=False)
    return dataframe

class voila_df(v.VuetifyTemplate):

    headers = List().tag(sync=True)
    items = List().tag(sync=True)
    template = Unicode().tag(sync=True)

    def __init__(self, df_):

        df = process_index(df_)

        super().__init__()
        columns = []
        v_slot = '<template v-slot:items="props">'
        for i, colname in enumerate(df.columns):
            if i == 0:
                list_item = {'text': colname, 'sortable': False, 'value': colname}
                v_slot = v_slot + f'<td>{{{{props.item.{colname}}}}}</td>'
            else:
                list_item = {'text': colname, 'value': colname, 'align': 'right'}
                v_slot = v_slot + f'<td class="text-xs-right">{{{{props.item.{colname}}}}}</td>'
            columns.append(list_item)
        v_slot = v_slot + '</template>'

        self.headers = columns
        self.items = df.to_dict(orient='records')
        self.template = '<v-data-table :headers="headers" :items="items" class="elevation-1">' \
            + v_slot \
            + '</v-data-table>'
