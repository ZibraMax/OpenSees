import opensees as ops
import math
import matplotlib.pyplot as plt

# --- Test OriHinge ---
ops.wipe()
ops.model('BasicBuilder', '-ndm', 3, '-ndf', 6)

# --- Parameters ---
L = 1.0
theta_0 = 210
phi = 2.0 * math.pi - theta_0 * math.pi / 180.0
x_coord = math.sin(math.pi / 3.0) * L
z_coord = math.sin(phi) * x_coord
x_coord2 = math.cos(phi) * x_coord

# --- Nodes ---
ops.node(1, 0.0, -0.5*L, 0.0)
ops.node(2, 0.0, 0.5*L, 0.0)
ops.node(3, x_coord, 0.0, 0.0)
ops.node(4, x_coord2, 0.0, z_coord)

# --- Supports ---
ops.fix(1, 1, 1, 1, 1, 1, 1)
ops.fix(2, 1, 1, 1, 1, 1, 1)
ops.fix(3, 1, 1, 1, 1, 1, 1)
ops.fix(4, 0, 0, 0, 1, 1, 1)

# --- Material ---
ops.uniaxialMaterial('Elastic', 1, 1000.0)

# --- Elements ---
ops.element('corotTruss', 1, 1, 2, 1.0, 1)
ops.element('corotTruss', 2, 2, 3, 1.0, 1)
ops.element('corotTruss', 3, 3, 1, 1.0, 1)
ops.element('corotTruss', 4, 1, 4, 1.0, 1)
ops.element('corotTruss', 5, 2, 4, 1.0, 1)
ops.element('OriHinge', 6, 3, 1, 2, 4, 0.3)

# --- Nodal Load ---
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
ops.load(4, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0)

# --- Static Analysis Setup ---
ops.system('Umfpack')
ops.constraints('Plain')
ops.numberer('RCM')
ops.test('NormDispIncr', 1.0e-8, 100)
ops.algorithm('Newton')
ops.integrator('ArcLength', 0.1, 1.0)
# ops.integrator('LoadControl', 0.1)
ops.analysis('Static')

# --- Perform Analysis ---
data = {}
data['step'] = []
data['load_factor'] = []
data['angle_degrees'] = []
data['disp_node4_z'] = []
print("Step,Load_Factor,Angle(degrees),Disp_Node4_Z")
M = 10
for i in range(1, M+1):
    if ops.analyze(1) != 0:
        print(f"Analysis failed at step {i}")
        break
    lam = ops.getLoadFactor(1)
    theta = ops.eleResponse(6, 'theta')
    disp_z = ops.nodeDisp(4, 3)
    print(f"{i},{lam},{theta[0] * 180.0 / math.pi},{disp_z}")
    data['step'].append(i)
    data['load_factor'].append(lam)
    data['angle_degrees'].append(theta[0] * 180.0 / math.pi)
    data['disp_node4_z'].append(disp_z)

plt.plot(data['angle_degrees'],data['load_factor'], 'o--')
plt.grid()
plt.xlabel('Angle (degrees)')
plt.ylabel('Load Factor')
plt.show()