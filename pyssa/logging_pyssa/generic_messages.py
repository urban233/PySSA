

def get_variables_values(filepath, step: str, variables_and_description: list[tuple[str, vars]]):
    messages = []
    for variable in variables_and_description:
        messages.append(f"Filepath: {filepath}; Step: {step}; {variable[0]}: {variable[1]}")
    return messages


def get_basic_variable_value(filepath, step, variable_name, variable):
    return f"Filepath: {filepath}; Step: {step}; {variable_name}: {variable}"
