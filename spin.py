
'''
import matplotlib.pyplot as plt

ax = plt.axes()
ax.arrow(0, 0, 0.2, 0.5, head_width=0.05, head_length=0.1,fc='k', ec='k')
plt.show()
'''

'''
import pylab as P
P.subplot(111)
# P.arrow( x, y, dx, dy, **kwargs )
P.arrow( 0.5, 0.8, 0.0, -0.2, fc="k", ec="k", fill = 'red',
head_width=0.05, head_length=0.1 )
P.show()
'''

import matplotlib.pyplot as plt
import matplotlib.patches as patches


f = open('configuration.txt')
content = []
for i in f:
    #print(i.rstrip().split())
    content.append(i.rstrip().split())
f.close()

print(content)

for i in range(1,len(content)):
    print(content[i][:-2])

T0 = patches.Arrow(
        0.4, 0.2, 0.0, 0.1, width=0.06,
        facecolor="red")

S0 = patches.Arrow(
        0.4, 0.2, 0.0, 0.1, width=0.06,
        facecolor="blue")

S11 = patches.Arrow(
        0.4, 0.2, 0.0, 0.1, width=0.06,
        facecolor="yellow")

S22 = patches.Arrow(
        0.4, 0.2, 0.0, 0.1, width=0.06,
        facecolor="green")

S1 = patches.Arrow(
        0.4, 0.2, 0.0, 0.1, width=0.06,
        facecolor="blue")

S2 = patches.Arrow(
        0.4, 0.2, 0.0, 0.1, width=0.06,
        facecolor="black")

arrows = []
for step in range(1,len(content)-2):
    for l,item in enumerate(content[step][:-2]):
        if step%2 == 1:
            if item == 'T0':
                T01 = patches.Arrow(
                0.2+l*0.15, 0.92-0.12*step+0.1, 0.0, 0.08, width=0.06,
                facecolor="red")
                arrows.append(T01)
                T02 = patches.Arrow(
                0.24+l*0.15, 1.0-0.12*step+0.1, 0.0, -0.08, width=0.06,
                facecolor="red")
                arrows.append(T02)
            elif item == 'S0':
                S01 = patches.Arrow(
                0.2+l*0.15, 0.92-0.12*step+0.1, 0.0, 0.08, width=0.06,
                facecolor="blue")
                arrows.append(S01)
                S02 = patches.Arrow(
                0.24+l*0.15, 1.0-0.12*step+0.1, 0.0, -0.08, width=0.06,
                facecolor="blue")
                arrows.append(S02)
            elif item == '11':
                Suu1 = patches.Arrow(
                0.2+l*0.15, 0.92-0.12*step+0.1, 0.0, 0.08, width=0.06,
                facecolor="yellow")
                arrows.append(Suu1)
                Suu2 = patches.Arrow(
                0.24+l*0.15, 0.9-0.12*step+0.1, 0.0, 0.08, width=0.06,
                facecolor="yellow")
                arrows.append(Suu2)
            elif item == '22':
                Sdd1 = patches.Arrow(
                0.2+l*0.15, 1.0-0.12*step+0.1, 0.0, -0.08, width=0.06,
                facecolor="black")
                arrows.append(Sdd1)
                Sdd2 = patches.Arrow(
                0.24+l*0.15, 1.0-0.12*step+0.1, 0.0, -0.08, width=0.06,
                facecolor="black")
                arrows.append(Sdd2)
            elif item == '1':
                S1 = patches.Arrow(
                0.2+l*0.15, 0.92-0.12*step+0.1, 0.0, 0.08, width=0.06,
                facecolor="green")
                arrows.append(S1)
            elif item == '2':
                S2 = patches.Arrow(
                0.2+l*0.15, 1.0-0.12*step+0.1, 0.0, -0.08, width=0.06,
                facecolor="brown")
                arrows.append(S2)
        else:
            if item == 'T0':
                T01 = patches.Arrow(
                0.15+l*0.15, 0.92-0.12*step+0.1, 0.0, 0.08, width=0.06,
                facecolor="red")
                arrows.append(T01)
                T02 = patches.Arrow(
                0.19+l*0.15, 1.0-0.12*step+0.1, 0.0, -0.08, width=0.06,
                facecolor="red")
                arrows.append(T02)
            elif item == 'S0':
                S01 = patches.Arrow(
                0.15+l*0.15, 0.92-0.12*step+0.1, 0.0, 0.08, width=0.06,
                facecolor="blue")
                arrows.append(S01)
                S02 = patches.Arrow(
                0.19+l*0.15, 1.0-0.12*step+0.1, 0.0, -0.08, width=0.06,
                facecolor="blue")
                arrows.append(S02)
            elif item == '11':
                Suu1 = patches.Arrow(
                0.15+l*0.15, 0.92-0.12*step+0.1, 0.0, 0.08, width=0.06,
                facecolor="yellow")
                arrows.append(Suu1)
                Suu2 = patches.Arrow(
                0.19+l*0.15, 0.9-0.12*step+0.1, 0.0, 0.08, width=0.06,
                facecolor="yellow")
                arrows.append(Suu2)
            elif item == '22':
                Sdd1 = patches.Arrow(
                0.15+l*0.15, 1.0-0.12*step+0.1, 0.0, -0.08, width=0.06,
                facecolor="black")
                arrows.append(Sdd1)
                Sdd2 = patches.Arrow(
                0.19+l*0.15, 1.0-0.12*step+0.1, 0.0, -0.08, width=0.06,
                facecolor="black")
                arrows.append(Sdd2)
            elif item == '1':
                if l==0:
                    S1 = patches.Arrow(
                    0.2+l*0.15, 0.92-0.12*step+0.1, 0.0, 0.08, width=0.06,
                    facecolor="green")
                else:
                    S1 = patches.Arrow(
                    0.12+l*0.15, 0.92-0.12*step+0.1, 0.0, 0.08, width=0.06,
                    facecolor="green")                    
                arrows.append(S1)
            elif item == '2':
                if l==0:
                    S2 = patches.Arrow(
                    0.2+l*0.15, 1.0-0.12*step+0.1, 0.0, -0.08, width=0.06,
                    facecolor="brown")
                else:
                    S2 = patches.Arrow(
                    0.12+l*0.15, 1.0-0.12*step+0.1, 0.0, -0.08, width=0.06,
                    facecolor="brown")                    
                arrows.append(S2)
            
            
    
fig6 = plt.figure()
ax6 = fig6.add_subplot(111, aspect='equal')
for p in arrows:
    ax6.add_patch(p)
for step in range(1,len(content)-2):
    print('energy:',content[step][-1])
    ax6.text(0.02, 1.05 - 0.12*step, str(content[step][-1])) 
    
'''
for p in [
    patches.Arrow(
        0.4, 0.2, 0.0, 0.1, width=0.06,
        facecolor="red"
    ),
    patches.Arrow(
        0.5, 0.2, 0.0, 0.1, width=0.06,
        facecolor="yellow"
    ),
]:
    ax6.add_patch(p)
'''
plt.show()
fig6.savefig('arrow6.png', dpi=90, bbox_inches='tight')



