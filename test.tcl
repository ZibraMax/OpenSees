# Test OriHinge
wipe
model BasicBuilder -ndm 3 -ndf 6

# --- Nodos ---
node 1 0.0 0.0 0.0
node 2 1.0 0.0 0.0
node 3 1.0 1.0 0.0
node 4 0.0 1.0 0.0

# --- Apoyo ---
fix 1 1 1 1 1 1
fix 2 1 1 1 1 1
fix 4 1 1 1 1 1

# --- Elemento OriHinge ---
# Sintaxis: element OriHinge tag nd1 nd2 nd3 nd4
element OriHinge 1 1 2 4 3

# --- Carga nodal ---
timeSeries Linear 1
pattern Plain 1 1 {
    load 3 0.0 0.0 10.0 0.0 0.0 0.0
}

# --- Análisis estático ---
system BandSPD
constraints Plain
numberer RCM
test NormDispIncr 1.0e-6 10
algorithm Newton
integrator LoadControl 1.0
analysis Static

analyze 1

# --- Imprimir desplazamientos ---
puts "Node 1: [nodeDisp 1]"
puts "Node 2: [nodeDisp 2]"
puts "Node 3: [nodeDisp 3]"
puts "Node 4: [nodeDisp 4]"
