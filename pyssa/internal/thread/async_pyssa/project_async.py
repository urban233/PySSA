

def search_for_not_matching_proteins(the_main_view_state, the_proteins) -> tuple:
    tmp_proteins = the_main_view_state.get_not_matching_proteins(the_proteins)
    return "result", tmp_proteins
