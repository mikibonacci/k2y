# Calcfunctions to perform the operations and store them in the AiiDA provenance/database.
import tempfile
import pathlib

from k2y.k2y import KcwQpDatabaseGenerator

from aiida.engine import calcfunction

from aiida import orm

@calcfunction
def generate_kcw_qp_database(
    yambo_retrieved,
    kcw_retrieved,
    QP_template_node=None
):
    """Generate a KcwQpDatabaseGenerator object and store it in the AiiDA database."""

    generator = KcwQpDatabaseGenerator.from_aiida(
        yambo_node_pk=yambo_retrieved.creator.pk,
        kcw_node_pk=kcw_retrieved.creator.pk,
        qp_template_node=QP_template_node,
    )
                
    generator.generate_mappings()
    # Create a temporary directory to store the database
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = pathlib.Path(temp_dir)
        new_db = generator.generate_QP_db_SinglefileData(filename="ndb.QP", temporary_dir=temp_path)
    
    return new_db