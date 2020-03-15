# coding: utf-8

'''Functions for generating a report via. papermill and nbconvert
'''

import logging
import nbformat
import os
import papermill as pm
import pkg_resources

from nbconvert import HTMLExporter

logger = logging.getLogger('root')

def generate_report(output_folder_path, nbargs):
    '''Generate an extraction report for a given run/output folder

    Parameters
    ----------
    output_folder_path: str
    nbargs: dict
        - Arguments to pass through to the notebook's parameters cell via. papermill

    Returns
    -------
    None

    '''

    report_template_path = pkg_resources.resource_filename('dart_id', 'notebooks/dart_id_report.ipynb')

    # Execute the notebook and write to disk
    notebook_output_path = os.path.join(output_folder_path, 'dart_id_report.ipynb')
    pm.execute_notebook(
        report_template_path,
        notebook_output_path,
        parameters=nbargs
    )

    # Convert notebook to HTML
    html_output_path = os.path.join(output_folder_path, 'dart_id_report.html')
    notebook_to_html(notebook_output_path, html_output_path)


def notebook_to_html(notebook_path, output_path):
    '''Convert a Jupyter/IPython notebook to an HTML file

    Parameters
    ----------
    notebook_path: str
    output_path: str

    Output
    ------
    None

    '''

    with open(notebook_path, 'r') as fp:
        notebook = nbformat.read(fp, nbformat.NO_CONVERT)
    
    # Instantiate the exporter
    html_exporter = HTMLExporter()
    # html_exporter.template_file = 'basic'

    # Exclude code cells from the HTML output
    html_exporter.exclude_input = True

    # Generate HTML
    (body, resources) = html_exporter.from_notebook_node(notebook)

    # print("Resources:", resources.keys())
    # print("Metadata:", resources['metadata'].keys())
    # print("Inlining:", resources['inlining'].keys())
    # print("Extension:", resources['output_extension'])

    with open(output_path, 'w') as fp:
        fp.write(body)
