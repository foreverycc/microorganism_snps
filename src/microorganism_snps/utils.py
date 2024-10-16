import logging

def conf_logger(log_file):
    '''
    Set up the logging system for the project.
    
    Parameters:
    - log_file: Path to the log file where logs will be saved
    
    Returns:
    - logger: Configured logger object
    '''
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    
    # File handler logs to the specified file
    lfh = logging.FileHandler(log_file)
    lfh.setFormatter(logging.Formatter('%(asctime)s:%(levelname)s:%(name)s: %(message)s'))
    logger.addHandler(lfh)
    
    # Stream handler logs to the console
    lsh = logging.StreamHandler()
    lsh.setFormatter(logging.Formatter('%(message)s'))
    logger.addHandler(lsh)
    
    return logger
