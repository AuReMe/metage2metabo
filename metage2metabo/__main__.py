from metage2metabo.m2m_workflow import run_workflow
import logging

logging.basicConfig(
    format='%(message)s', level=logging.DEBUG)  #%(name)s:%(message)s
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

def main_workflow(args=None):
    run_workflow()

if __name__ == "__main__":
    main_workflow()