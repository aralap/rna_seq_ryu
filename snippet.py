def get_config_parameter(parameter):
    config_file = open('config.cfg','r')
    for line in config_file:
        if not line.startswith('#') and line != '\n':
            line = line.split()
            if parameter in line[0]:
                config_file.close()
                return line[2].strip('\'')