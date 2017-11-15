function print_mess     
comps = graphconncomp(handles.SimGraph, 'Directed', false);
    msg = [msg sprintf(' - %d connected components found', comps)];