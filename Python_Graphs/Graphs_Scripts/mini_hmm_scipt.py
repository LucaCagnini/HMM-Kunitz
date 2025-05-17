from graphviz import Digraph

# Create graphs
dot = Digraph(format='png')
dot.attr(rankdir='TR', size='10,6')
dot.node("Start", shape="circle", style="filled", fillcolor="lightgray")
dot.node("End", shape="circle", style="filled", fillcolor="lightgray")

# example states for three columns
states_M = ['M1', 'M2', 'M3','M4']
states_I = ['D1', 'D2', 'D3','D4']
states_D = ['I1', 'I2', 'I3','I4']

# nodes
for i in range(0, 4):
    dot.node(states_M[i], shape='square', style="filled", fillcolor="lightblue")
    dot.node(states_I[i], shape="circle", style="filled", fillcolor="lightgreen")
    dot.node(states_D[i], shape="diamond", style="filled", fillcolor="orange")




with dot.subgraph() as s:
    s.attr(rank='same')
    s.node('M1')
    s.node('M2')
    s.node('M3')
    s.node('M4')
    #s.node("Start")
    #s.node("End")

with dot.subgraph() as s:
    s.attr(rank='same')
    s.node('D1')
    s.node('D2')
    s.node('D3')
    s.node('D4')

with dot.subgraph() as s:
    s.attr(rank='same')
    s.node('I1')
    s.node('I2')
    s.node('I3')
    s.node('I4')

# Start transitions
dot.edge('Start', 'I1', label="Start")
dot.edge('Start', 'M1', label="Start")
dot.edge('Start', 'D1', label="Start")



#M1
dot.edge('M1', 'M2', label='m1->m2')
dot.edge('M1', 'I1', label='m1->i1')
dot.edge("M1", "D2", label="m1 →d2")
#I1
dot.edge('I1', 'I1', label='i1->i1')
dot.edge('I1', 'D2', label='i1->d1')
dot.edge('I1', 'M2', label='i1->m2')
#D1
dot.edge('D1', 'D2', label='d1->d2')
dot.edge('D1', 'I1', label='d1->i1')
dot.edge('D1', 'M2', label='d1->i1')


#M2
dot.edge('M2', 'M3', label='m2->m3')
dot.edge('M2', 'I2', label='m2->i2')
dot.edge('M2', 'D3', label='m2->d3')
#I2
dot.edge("I2", "D3", label="i2->d3")
dot.edge("I2", "M3", label="i2 → i3")
dot.edge("I2", "I2", label="i2 → i2")
#D2
dot.edge("D2", "D3", label="d2 → d3")
dot.edge("D2", "M3", label="m2 → m3")
dot.edge('D2', 'I2', label='d2 →i2')

#M3
dot.edge('M3', 'M4', label='m3->m4')
dot.edge('M3', 'I3', label='m3->i3')
dot.edge('M3', 'D4', label='m3->d4')
#I3
dot.edge("I3", "D4", label="i3->d4")
dot.edge("I3", "M4", label="i3 → i4")
dot.edge("I3", "I3", label="i3 → i3")
#D2
dot.edge("D3", "D4", label="d3 → d4")
dot.edge("D3", "M4", label="m3 → m4")
dot.edge('D3', 'I3', label='d3 →i3')


#M4
dot.edge("M4", "I4", label="m4 → i4")
dot.edge("M4","End", label="m4 → End")
#D3
dot.edge("D4", "End", label="d4 → End")
#I3
dot.edge("I4", "I4", label="i4 → i4")
dot.edge("I4", "End", label="i4 → End")


# Save
dot.render('mini_hmm')
