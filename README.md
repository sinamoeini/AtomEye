
# AtomEye: atomistic configuration viewer

<p align="center">


[Gallery](gallery.html) | [Features](#features) | [File formats](#formats) | [Manual](#manual) | [Mac Keybinding](#MacKeybinding) | [Download](#download) | [Utilities](utils.html) | [FAQ](#FAQ) | [Bug report](#bugreport) | [History](#history)

</p>
<font color="red">Heavy users should try the parallelized [AtomEye version 3 (2012)](../A3/A3.html), with scripting capability.</font>

AtomEye will always be free. Please [cite](http://scholar.google.com/scholar?cites=6489407935302268968)
[J. Li, _Modelling Simul. Mater. Sci. Eng._ **11** (2003) 173](Doc/Li03a.pdf)
if you use figures and movies produced by AtomEye, so your colleagues may know about this free tool. Thanks.



* * *

<a name="features">

## Features

*   order-N in both execution time and memory used, where N is the number of atoms; designed for condensed-matter systems, no problem with >1 million atoms
*   auto-detect 8, 16 and 32-bit shared memory or remote X-displays
*   multiple resizable windows in POSIX threads; 0% CPU usage if not moving
*   geometrically exact area-weighted antialiasing for atoms, bonds, and wireframes
*   fast rendering of atoms by caching pixel- and z-maps in the main memory
*   quick toggle between parallel and perspective projections
*   full 3-D navigation
*   support periodic boundary conditions
*   support PDB input file format
*   support arbitrary-precision and extendable CFG input file format for large-scale, reloadable molecular dynamics simulations
*   auto-decompress gzip- or bzip2-compressed input configuration files
*   JPEG, PNG and high-resolution EPS screenshots
*   customizable atom radii and coloring schemes
*   coordination number color-encoding with customizable cutoff radii and invisibility controls
*   local atomic von Mises shear strain invariant color-encoding
*   user-defined property color-encoding, in hsv, jet, and other colormapping choices
*   color-marking an initial configuration to track subsequent atomic displacements
*   cooperative X-terminal input with GNU readline / history
*   up to 16 arbitrary cutting planes with advancing / rotation / flipping controls
*   animation script for making movies

* * *

</a><a name="formats">

## Configuration file formats

AtomEye works with atomistic configurations, the format of which should follow the general guidelines below from our experience:

</a>

<a name="formats"></a>
*   <a name="formats">One file stores only one configuration.</a> [XMol](http://neon.orch.ruhr-uni-bochum.de/progs/xmol.html) XYZ's concatenated format is not recommended.
*   For portability, the file should be in ASCII plain text that a human can read.
*   There should _minimally_ be chemical symbol designation for each atom and some representation for its position.
*   It is recommended to define comment line in the format, and a configuration file should be self-explanatory with the aid of comments.

AtomEye currently supports the following formats:   [PDB](#PDB) | [standard CFG](#standard_CFG) | [extended CFG](#extended_CFG)<a name="compressed_config_file"></a>

<a name="compressed_config_file">All configuration files can be compressed by</a> [gzip](http://www.gzip.org/) / [bzip2](http://sources.redhat.com/bzip2/) and then fed directly into AtomEye, which looks for magic number at the first few bytes, and if found would automatically decompress it by calling shell gzip/bzip2\. In order to use this feature, the user must install gzip/bzip2 executable in his shell PATH.

[Standard](#standard_CFG) or [extended](#extended_CFG) CFG formats are _strongly_ recommended over the [PDB](#PDB) format for hassle-free rendering by AtomEye. Among other things, [standard](#standard_CFG) or [extended](#extended_CFG) CFG formats allow arbitrary numerical precision. Also, enforcing PBC is quite a bit of work for the [PDB](#PDB) format, but is built-in in [standard](#standard_CFG) or [extended](#extended_CFG) CFG formats.

<a name="PDB"></a>

#### [Protein Data Bank](http://www.rcsb.org/pdb/) format

This format is widely used for storing atomic structure. Its specifications can be found [here](http://www.rcsb.org/pdb/docs/format/pdbguide2.2/guide2.2_frame.html). The standard renderer of PDB file is [Rasmol](http://www.bernstein-plus-sons.com/software/rasmol/).

PDB examples:    [DNA.pdb](Gallery/DNA/DNA.pdb) | [Nanotube8x3x1.pdb](Gallery/Nanotube/Nanotube8x3x1.pdb) | [SiShuffle.pdb](Gallery/SiShuffle/SiShuffle.pdb)

In order to enforce periodic boundary condition (PBC) behavior in AtomEye, one must insert the following two lines at the beginning of the file:

<pre>HEADER  500 MUST_PBC
CRYST1   13.303   37.627    7.681  60.00  70.00  80.00 P 1           1
</pre>

in which the first line is a secret code that stays; the second line contains six numbers,

*   _a=13.303_:  the length of the first box edge in Angstrom
*   _b=37.627_:  the length of the second box edge in Angstrom
*   _c=7.681_:  the length of the third box edge in Angstrom
*   _alpha=60.00_:  the angle between the second and third box edges in degrees
*   _beta=70.00_:  the angle between the first and third box edges in degrees
*   _gamma=80.00_:  the angle between the first and second box edges in degrees

which should be modified accordingly. The box edges are taken to be _h1=(u1,0,0)_, _h2=(u2,u3,0)_, _h3=(u4,u5,u6)_ in Cartesian coordinates, with _u1,u2,u3,u4,u5,u6_ having one to one correspondence with _a,b,c,alpha,beta,gamma_. This connection is important because the atom positions will be specified in Cartesian coordinates, so one must unequivocally specify the three PBC edges in Cartesian coordinates too, fixing the rotational degrees of freedom.

A configuration file is considered to be in PDB format if there is '.pdb' or '.PDB' in the filename.<a name="standard_CFG"></a>

#### <a name="standard_CFG">Standard CFG format</a>

<a name="standard_CFG">

The standard CFG format ensures seamless transition and _complete_ information passage from one MD simulation to the other. Thus particle velocities _must_ be specified, and floating-point numbers are usually saved to the 16th significant digit. Aside from an atom's position and velocity, properties such as the local energy are deemed potential-dependent and therefore cannot be specified in the standard format.

The shortcoming of the standard CFG format is its lack of extensibility, and the filesize tends to be noticeably larger for very large configurations (>200K atoms).

</a>

<a name="standard_CFG">Standard CFG examples:</a> [SiVacancy.cfg](Gallery/SiVacancy/SiVacancy.cfg) | [Si_screw_dipole.cfg](Gallery/Si_screw_dipole/Si_screw_dipole.cfg) | [Cu_NanoXtal.cfg.bz2](Gallery/Cu_NanoXtal/Cu_NanoXtal.cfg.bz2)

Take a look at [SiVacancy.cfg](Gallery/SiVacancy/SiVacancy.cfg). The format consists of system specifications and atom specifications. The system specifications are of the form 'tagname = value units', where 'tagname' and 'units' are fixed, and 'value' needs to be filled in. Some tagnames are _required_, such as 'Number of particles', 'H0(1,1)' to 'H0(3,3)'. Some are optional, such as 'A', 'eta(1,1)', 'R'. A line starting with '#' is comment. Here the comments explain the tags immediately above, which is a good habit.<a name="standard_CFG_example"></a>

<pre><a name="standard_CFG_example">Number of particles = 1727
# (required) this must be the first line</a> </pre>

<a name="standard_CFG_example">

<pre>A = 1.0 Angstrom (basic length-scale)
# (optional) basic length-scale: default A = 1.0 [Angstrom]
</pre>

<pre>H0(1,1) = 32.5856986704313 A
H0(1,2) = 0 A
H0(1,3) = 0 A
# (required) this is the supercell's 1st edge, in A

H0(2,1) = 8.64689152483509e-16 A
H0(2,2) = 32.5856986704313 A
H0(2,3) = 0 A
# (required) this is the supercell's 2nd edge, in A

H0(3,1) = 8.64689152483509e-16 A
H0(3,2) = 8.64689152483509e-16 A
H0(3,3) = 32.5856986704313 A
# (required) this is the supercell's 3rd edge, in A
</pre>

<pre>Transform(1,1) = 1
Transform(1,2) = 0
Transform(1,3) = 0
Transform(2,1) = 0
Transform(2,2) = 1
Transform(2,3) = 0
Transform(3,1) = 0
Transform(3,2) = 0
Transform(3,3) = 1
# (optional) apply additional transformation on H0:  H = H0 * Transform;
# default = Identity matrix.

eta(1,1) = 0
eta(1,2) = 0
eta(1,3) = 0
eta(2,2) = 0
eta(2,3) = 0
eta(3,3) = 0
# (optional) apply additional Lagrangian strain on H0:
# H = H0 * sqrt(Identity_matrix + 2 * eta);
# default = zero matrix.
</pre>

<pre># ENSUING ARE THE ATOMS, EACH ATOM DESCRIBED BY A ROW
# 1st entry is atomic mass in a.m.u.
# 2nd entry is the chemical symbol (max 2 chars)

# 3rd entry is reduced coordinate s1 (dimensionless)
# 4th entry is reduced coordinate s2 (dimensionless)
# 5th entry is reduced coordinate s3 (dimensionless)
# real coordinates x = s * H,  x, s are 1x3 row vectors

# 6th entry is d(s1)/dt in basic rate-scale R
# 7th entry is d(s2)/dt in basic rate-scale R
# 8th entry is d(s3)/dt in basic rate-scale R
R = 1.0 [ns^-1]
# (optional) basic rate-scale: default R = 1.0 [ns^-1]
</pre>

<pre>28.0855 Si .0208333333333333 .0208333333333333 .0208333333333333 0 0 0
28.0855 Si .0625 .0625 .0625 0 0 0
28.0855 Si .0208333333333333 .104166666666667 .104166666666667 0 0 0
28.0855 Si .0625 .145833333333333 .145833333333333 0 0 0
28.0855 Si .104166666666667 .0208333333333333 .104166666666667 0 0 0
28.0855 Si .145833333333333 .0625 .145833333333333 0 0 0
....
</pre>

The atom specifications are the rows above, one row for each atom, with the meaning of each column entry explained in the comments. The reduced coordinates s1,s2,s3 should be between 0 and 1\. The Cartesian coordinates of the atoms can be determined from the reduced coordinates as,

<pre> x  =  s1 * H(1,1)  +  s2 * H(2,1)  +  s3 * H(3,1)
y  =  s1 * H(1,2)  +  s2 * H(2,2)  +  s3 * H(3,2)
z  =  s1 * H(1,3)  +  s2 * H(2,3)  +  s3 * H(3,3)
</pre>

</a>

<a name="standard_CFG_example">'A', 'R' are the lengthscale and ratescale of H0[][] and ds[]/dt, respectively, which by default take the values Angstrom and nanosecond^-1\. They can be altered to manually dilate or heat up the system. 'eta[][]' is the optional Lagrangian strain which is 0 by default, using which we can apply an additional deformation on H0[][] to get the actual box shape H[][]. These are amenities that entry-level users can do without.</a><a name="extended_CFG"></a>

#### <a name="extended_CFG">Extended CFG format</a>

<a name="extended_CFG"></a>

<a name="extended_CFG">Extended CFG examples:</a> [BubbleRaftBefore.cfg](Gallery/BubbleRaft/BubbleRaftBefore.cfg) | [BubbleRaftAfter1.cfg](Gallery/BubbleRaft/BubbleRaftAfter1.cfg)

Extended CFG format addresses the problem of filesize and extensibility. Velocity data is no longer required for the atoms, and floating-point precision can be what one sees fit for the application instead of the recommended 16 significant digits in case of the standard CFG format. Also, atomic mass and chemical symbol are now assigned on a block-by-block basis instead of on an atom-by-atom basis. All standard CFG system tags are still valid in the extended CFG format, with addition of the 'entry_count' keyword, whose appearance determines whether this is standard or extended CFG format.

<pre>Number of particles = 16200
A = 4.37576470588235 Angstrom (basic length-scale)
H0(1,1) = 127.5 A
H0(1,2) = 0 A
H0(1,3) = 0 A
H0(2,1) = 0 A
H0(2,2) = 119.501132067411 A
H0(2,3) = 0 A
H0(3,1) = 0 A
H0(3,2) = 0 A
H0(3,3) = 3 A
.NO_VELOCITY.
entry_count = 11
auxiliary[0] = kine [reduced unit]
auxiliary[1] = pote [reduced unit]
auxiliary[2] = s11 [reduced unit]
auxiliary[3] = s22 [reduced unit]
auxiliary[4] = s12 [reduced unit]
auxiliary[5] = hydro [reduced unit]
auxiliary[6] = mises [reduced unit]
auxiliary[7] = Lmin [reduced unit]
1.000000
Ar
0.005 0.01232 0.5 0 -2.9819 0.79705 2.3326 0.022255 1.5648 0.76808 102.82
....
</pre>

In above, 'entry_count = 11' means there are 11 entries per row for an atom. To omit velocity data in rows, write '.NO_VELOCITY.' line before the 'entry_count = 11' line. In this case because 11=3+8, we must have 8 so-called _auxiliary properties_ per atom. If on the other hand '.NO_VELOCITY.' does not appear before 'entry_count = 11', then the velocities (ds[]/dt in R) will be specified after s[] in each row, and since 11=3+3+5, we will have only 5 auxiliary properties.

For AtomEye to know what these auxiliary properties are, so it can provide the correct help information, one should also fill in 'auxiliary[i] = name [unitname]' lines, where 'i' runs from 0 to 8-1=7 here, 'name' is a single word for the property's name, and unitname is the name of the unit in which the values are given. For example, 'auxiliary[3] = mises [GPa]' would mean that the _fourth_ auxiliary property, which is the 7th column (or the 10th column if '.NO_VELOCITY.' clause is not given), are the von Mises stress invariants in GPa.

Lastly, atomic mass and chemical symbol are now assigned on a block-by-block basis. In the above file _all_ atoms would be 1.000000 amu Ar, which is a pseudonym for soap bubble. Apparently this format can save considerable space for monoatomic configuration. When we have binary compounds, we can write,

<pre>28.0855
Si
0.005 0.01232 0.5 0 -2.9819 0.79705 2.3326 0.022255 1.5648 0.76808 102.82
12.011
C
0.091667 0.01232 0.5 0 -2.9491 1.3535 3.9801 0.9655 2.6668 1.63 103.07
28.0855
Si
0.31833 0.01232 0.5 0 -2.9787 0.84615 2.4817 0.45222 1.6639 0.93449 103.02
12.011
C
0.765 0.01232 0.5 0 -2.9431 1.4263 4.1825 -1.0568 2.8044 1.7367 102.97
</pre>

which takes up no more space than the standard CFG format under UNIX. Of course a space-saving way is to write the above as,

<pre>28.0855
Si
0.005 0.01232 0.5 0 -2.9819 0.79705 2.3326 0.022255 1.5648 0.76808 102.82
0.31833 0.01232 0.5 0 -2.9787 0.84615 2.4817 0.45222 1.6639 0.93449 103.02
12.011
C
0.091667 0.01232 0.5 0 -2.9491 1.3535 3.9801 0.9655 2.6668 1.63 103.07
0.765 0.01232 0.5 0 -2.9431 1.4263 4.1825 -1.0568 2.8044 1.7367 102.97
</pre>

* * *

Sometime it is more convenient to have the auxiliary properties in a stand-alone '.aux' file, which can be patched to the current application on-demand (press 'F11'). This file should look like,

<pre>102.54352   3.781
-54.324521 -9.035
22.870594   0.785
-9.8543543  6.834
....
</pre>

where each line contains the property values for one atom, and the total number of lines is equal to the total number of atoms. An example would be [amorphous_Si.cfg](Gallery/Amorphous_Si/amorphous_Si.cfg)  |  [amorphous_Si.aux](Gallery/Amorphous_Si/amorphous_Si.aux) file combination. See [auxiliary property coloring](#auxiliary_property_coloring) for how the coloring intensity is controlled.

* * *

<a name="manual">

## Manual

</a>

[usage](#usage) | [NumLock](#NumLock) | [help key](#help_key) | [rotate object](#rotate_object) | [anchor control](#anchor_control) | [toggle bond mode](#toggle_bond_mode) | [toggle parallel / perspective projection](#toggle_perspective) | [scale atom radii](#scale_atom_radii) | [change bond radius](#change_bond_radius) | [toggle wireframe mode](#toggle_wireframe_mode) | [cutoff control](#cutoff_control) | [upright viewframe](#upright_viewframe) | [inquire atom info](#inquire_atom_info) | [inquire geometrical info](#inquire_geometrical_info) | [pull viewport closer/away from anchor](#pull_closer_away) | [shift object](#shift_object) | [changing gear](#changing_gear) | [coordination number coloring](#coordination_number_coloring) | [shear strain coloring](#shear_strain_coloring) | [central symmetry coloring](#central_symm_coloring) | [least-square atomic local strain coloring](#least_square_strain_coloring) | [auxiliary property coloring](#auxiliary_property_coloring) | [extra patch file coloring](#extra_color_patch) | [make atoms invisible](#make_atoms_invisible) | [make bond invisible](#make_bond_invisible) | [change atom color and/or radius](#change_atom_color_radius) | [change bond color and/or radius](#change_bond_color_radius) | [shift object under PBC](#shift_object_under_PBC) | [change background color](#change_bgcolor) | [change view angle amplification](#change_total_view_angle) | [new / clone / quit viewport](#new_clone_quit_viewport) | [print system status](#system_status) | [find atom](#find_atom) | [goto position](#goto_position) | [resize window](#resize_window) | [jpg screenshot](#making_jpg) | [png screenshot](#making_png) | [eps screenshot](#making_eps) | [toggle auto shell viewer](#toggle_viewer) | [load new config](#load_new_config) | [reload config](#reload_config) | [sequential config list browsing](#step_config) | [define and trace color tiling blocks](#color_tiling) | [creating and manipulating cutting planes](#cutting_planes) | [save atom indices in a file](#save_atom_indices) | [making movie](#making_movie)

<a name="usage">

#### Usage:

</a>[Download](#download) binary from browser and save as A

*   % chmod 755 A
*   % ./A filename

You must ensure that the "xterm" command is in your shell PATH and can be directly called. See also [here](#compressed_config_file) if your configuration file has been compressed.<a name="NumLock">#### NumLock:

In general, your NumLock key should be in the 'off' state. Otherwise, AtomEye will regard it as equivalent to having CapsLock 'on', which is equivalent to having</a> [Meta](#Meta)+ modifier for every key you press.<a name="help_key">

#### Help Key:

Press 'F1' or 'h' for help.</a><a name="rotate_object">

#### Rotate Object (or so you think):

</a>

<a name="rotate_object">This can be accomplished either with arrow keys or mouse drag. Press Left, Right, Up, Down keys to rotate object as if you are rolling a crystal ball with the object embedded in it. For in-plane rotation, use Shift+Up (clockwise) and Shift+Down (counter-clockwise). The angle of rotation corresponding to one keystroke is controlled by the</a> [gearbox value](#changing_gear) multiplied by pi. Therefore, at gear-9 corresponding to gearbox value of 0.5, each keystroke will flip the object by pi/2, which comes in handy when viewing large configurations.

One can also rotate the viewport by pressing left button and drag the mouse. The "crystal ball" can be imagined as being hinged at the viewport center and whose diameter is somewhat smaller than the viewport. If the pointer should fall onto the crystal ball's surface, the rotation follows a geodesic path connecting A and B on the ball's surface. On the other hand, if the pointer falls outside of the ball, then dragging the mouse will just rotate the object clockwise/counter-clockwise in-plane.

It can be foreseen that when the configuration is large, it is much more efficient to rotate the viewport than to rotate the object. In fact, contrary to what the heading suggests, this _is_ how I implement the rotation. To get the correct look and feel, one must also _translate_ the viewport after the rotation. It is done by the following rule: we stipulate that there exists a location in space whereby the combined effect of rotation and translation would not change its _subjective_ position with respect to the viewport. This location is called the _anchor_, named so because when we "rotate the object", it appears to be hinging on the anchor. Please visit [anchor control](#anchor_control) to see how to change the anchor position.<a name="anchor_control"></a>

#### <a name="anchor_control">Anchor Control:</a>

<a name="anchor_control">

At the start of a session, the anchor is taken to be the center of the box, so the object appears to be hinged at the box center if rotated. This state can be recovered whenever 'w' is pressed.

On the other hand, right-clicking on an atom transfers the anchor to that atom's position. This allows a streamlined right click + drag action that pulls the viewport closer to any atom you want to see in detail. Also, once you right-clicked on an atom, ensuing rotations will be hinged on that atom, which is convenient for studying local atomic arrangements.

</a>

<a name="anchor_control">If you are in the</a> [bond mode](#toggle_bond_mode), then right-clicking on a bond will also set the anchor position to the bond center.

Press 'a' to shift the viewport so the anchor is seen right at the middle of the viewport.<a name="toggle_bond_mode"></a>

#### <a name="toggle_bond_mode">Toggle Bond Mode:</a>

<a name="toggle_bond_mode"></a>

<a name="toggle_bond_mode">Press 'b' to toggle whether to draw bonds or not. To change the bonding cutoff between two species of atoms, visit</a> [Cutoff Control](#cutoff_control).<a name="toggle_perspective"></a>

#### <a name="toggle_perspective">Toggle Between Parallel / Perspective Projection:</a>

<a name="toggle_perspective"></a>

<a name="toggle_perspective">Press 'Tab' to toggle between parallel and perspective projection rendering methods. Parallel projection is a limiting case of perspective projection, where the viewpoint is very far away from the object but the</a> [view angle](#change_total_view_angle) is also turned very small. It is useful for discerning atomic rows and planes.<a name="scale_atom_radii"></a>

#### <a name="scale_atom_radii">Scale Atom Radii:</a>

<a name="scale_atom_radii"></a>

<a name="scale_atom_radii">Press '</a>[PageUp](#PageUp)' to increase atom radii and '[PageDown](#PageDown)' to decrease atom radii rendered on screen by a common factor. The rate of change is controlled by the [gearbox](#changing_gear).<a name="change_bond_radius"></a>

#### <a name="change_bond_radius">Change Bond Radius:</a>

<a name="change_bond_radius"></a>

<a name="change_bond_radius">Press 'Home' to increase bond radius and 'End' to decrease bond radius drawn on the screen. The rate of change is controlled by the</a> [gearbox](#changing_gear).<a name="toggle_wireframe_mode"></a>

#### <a name="toggle_wireframe_mode">Toggle Wireframe Mode:</a>

<a name="toggle_wireframe_mode"></a>

<a name="toggle_wireframe_mode">Press 'i' to toggle wireframe mode that determines how to render the PBC box. You have the freedom to select no wireframe, monochromatic wireframe, RGB wireframe, and so on. In the case of the RGB wireframe, the three axes correspond to (s1,0,0), (0,s2,0), (0,0,s3), respectively; the point they meet corresponds to the origin (0,0,0).</a><a name="cutoff_control"></a>

#### <a name="cutoff_control">Cutoff Control:</a>

<a name="cutoff_control">

Press 'r' for cutoff control, which decides how close two atoms need to be with each other to be considered nearest neighbor (if they are, then each atom's coordination number will be increased by 1, and a bond will be drawn between them if under the bond mode).

</a>

<a name="cutoff_control">You will be first inquired about which two species. For example, if you want to change the cutoff distance between silicon and carbon, then you should enter "Si C" and press enter. Then, with each Ctrl+Home/Ctrl+End, the cutoff distance will be increased/decreased; the rate of change is controlled by the</a> [gearbox](#changing_gear). When you are satisfied, you should press 'r' again to close the change, so you can't modify two pairing types at the same time.<a name="upright_viewframe"></a>

#### <a name="upright_viewframe">Make Viewframe Upright:</a>

<a name="upright_viewframe"></a>

<a name="upright_viewframe">Sometimes too much rotation is a bad thing, and you want the viewframe upright again just like at the beginning. Press 'u'.</a><a name="inquire_atom_info"></a>

#### <a name="inquire_atom_info">Inquire Atom Information:</a>

<a name="inquire_atom_info"></a>

<a name="inquire_atom_info">Right click on an atom to display relevant information about it in the xterm. This will also set the</a> [anchor](#anchor_control) to this atom.

Right-clicking on a bond will display information about this bond. This will set the [anchor](#anchor_control) to the bond center.<a name="inquire_geometrical_info"></a>

#### <a name="inquire_geometrical_info">Inquire Geometrical Information:</a>

<a name="inquire_geometrical_info"></a>

<a name="inquire_geometrical_info">AtomEye remembers the last four atoms that you have</a> [clicked on](#inquire_atom_info).

*   Pressing ',' (comma) would print out the Cartesian separation and distance between the last two atoms clicked.
*   Pressing '.' (period) would print out the bond angle between the last three atoms clicked.
*   Pressing '/' (slash) would print out the dihedral angle between the last four atoms clicked.

<a name="pull_closer_away">

#### Pull Viewport Closer/Away from the Anchor:

</a>

<a name="pull_closer_away">Right click in the window and drag mouse to pull closer. Alternatively, use IMWheel (or Ctrl+IMWheel for quicker action) to pull viewport closer/away from the anchor. If you click on an atom, that atom will automatically become the</a> [anchor](#anchor_control). If you click on a bond, the bond center will also become the [anchor](#anchor_control).<a name="shift_object"></a>

#### <a name="shift_object">Shift Object (or so you think):</a>

<a name="shift_object"></a>

<a name="shift_object">Ctrl+Left, Ctrl+Right, Ctrl+Up, Ctrl+Down will shift the object in plane. Ctrl+Shift+Up will send the object further from viewport, Ctrl+Shift+Down will pull it closer.</a><a name="changing_gear"></a>

#### <a name="changing_gear">Changing Gearbox Value:</a>

<a name="changing_gear">

Press numeral keys '0' to '9', and the gearbox value will be switched to,

<pre>[1]0.001 [2]0.002 [3]0.005 [4]0.010 [5]0.020
[6]0.050 [7]0.100 [8]0.200 [9]0.500 [0]0.150
</pre>

The gearbox value controls all rate of change, such as angle of rotation, amount of translation, the rate of atom radius change, etc.</a><a name="coordination_number_coloring">

#### Coordination Number Coloring:

</a>

<a name="coordination_number_coloring">Press 'k' to toggle coordination number coloring. Coordination number is an empirical measure of how many nearest neighbors (could be of various species) there are for a particular atom. The definition of "nearest neighbors" can be changed by</a> [cutoff control](#cutoff_control).

To clearly see defect cores, you often need to remove the perfectly coordinated atoms. This is done in [Make Atoms Invisible](#make_atoms_invisible).<a name="shear_strain_coloring"></a>

#### <a name="shear_strain_coloring">Atomistic Local von Mises Shear Strain Invariant Coloring:</a>

<a name="shear_strain_coloring"></a>

<a name="shear_strain_coloring"></a>[Meta](#Meta)+g will color-encode the atoms according to their [local von Mises shear strain invariant](Doc/vonMisesInvariant.pdf). Shift+g will toggle the flag controlling whether to subtract off the system-averaged strain tensor or not before computing the invariant; the default is no. The controls of colormap, visibilities etc. are identical to that of [auxiliary properties coloring](#auxiliary_property_coloring).<a name="central_symm_coloring"></a>

#### <a name="central_symm_coloring">Central Symmetry Parameter Coloring:</a>

<a name="central_symm_coloring"></a>

<a name="central_symm_coloring"></a>[Meta](#Meta)+h will color-encode the atoms according to their [central symmetry parameter](Doc/CentralSymmetry.pdf) _c_'s. Shift+h will prompt the user to change the maximum number of neighbors _M_ used in the computation; the default being the most popular coordination number of the configuration rounded even. The controls of colormap, visibilities etc. are identical to that of [auxiliary properties coloring](#auxiliary_property_coloring).

Because calculating the central symmetry parameters requires a [neighborlist](Doc/neighborlist.pdf) without pairwise saving, it is not calculated by default when the configuration is loaded, but only computed after [Meta](#Meta)+h is pressed. _c_ is between [0,1], and its average is usually less than 0.5 even for amorphous structures. An intrinsic stacking fault in FCC crystal would possess two layers of atoms with _c_'s at about 0.042, and a perfect crystal should have _c_'s less than 0.01 even with thermal fluctuations. So a recommended threshold value for visualizing planar faults is 0.01, which can be set by shifting the [gearbox](#changing_gear) to level-4 and pressing Ctrl+[PageUp](#PageUp) once. Also, as [auxiliary properties coloring](#auxiliary_property_coloring) explains, one needs to press 'Ctrl+T' to keep the above set thresholds for ensuing configurations. See also [auxiliary properties coloring](#auxiliary_property_coloring) for further controls.<a name="least_square_strain_coloring"></a>

#### <a name="least_square_strain_coloring">Least-Square Atomic Local Strain Tensor Coloring:</a>

<a name="least_square_strain_coloring">We have implemented</a> [atomic local strain tensor](annotate_atomic_strain/Doc/main.pdf) coloring as a core functionality of AtomEye (please cite F. Shimizu, S. Ogata and J. Li, "[Theory of Shear Banding in Metallic Glasses and Molecular Dynamics Calculations](../../Papers/07/Shimizu07a.pdf)," _Materials Transactions_ **48** (2007) 2923-2927, if you use this characterization). It is different from the [von Mises shear strain invariant coloring](#shear_strain_coloring) in that two configurations are required, one current, and one reference configuration. It is also much more powerful and robust than the former because it does not rely on assumptions of prior high lattice symmetry.

To use, you need to load in a configuration first, and then press 'Esc' key. This would _imprint_ the present coordinates as the reference. Any configuration you [load](#load_new_config) or [step](#step_config) afterward would automatically trigger local strain calculations, the results of which are appended as [auxiliary properties](#auxiliary_property_coloring). Right-click on an atom to see the full local transformation matrix J, defined as

<pre> dx_{ij}(now)  \approx   dx_{ij}(reference) J
</pre>

where dx's are row vectors, and j are i's [nearest neighbors](#cutoff_control) (in the present configuration, as we don't want to store two neighborlists).

Note that for the results to be meaningful, the two configuration files must be _isoatomic_, meaning they must have the same set of atoms, with the same indexing (order). All that that may change between the two configurations are the [H[][] matrix and s[] coordinates](#standard_CFG).

This functionality is also spun off as a standalone utility called [annotate_atomic_strain](utils.html#annotate_atomic_strain). Please cite F. Shimizu, S. Ogata and J. Li, "[Theory of Shear Banding in Metallic Glasses and Molecular Dynamics Calculations](../../Papers/07/Shimizu07a.pdf)," _Materials Transactions_ **48** (2007) 2923-2927, if you use this utility.<a name="auxiliary_property_coloring"></a>

#### <a name="auxiliary_property_coloring">Auxiliary Property Coloring:</a>

<a name="auxiliary_property_coloring"></a>

<a name="auxiliary_property_coloring"></a>[Meta](#Meta)+[0-9,a-f] will color-encode [auxiliary properties](#extended_CFG) 0 to 15, whereas [Meta](#Meta)+[0-9,a-f] with CapsLock ON will color-encode [auxiliary properties](#extended_CFG) 16 to 31\. See [BubbleRaftBefore.cfg](Gallery/BubbleRaft/BubbleRaftBefore.cfg), [BubbleRaftBefore.jpg](Gallery/BubbleRaft/BubbleRaftBefore.jpg), [BubbleRaftAfter1.cfg](Gallery/BubbleRaft/BubbleRaftAfter1.cfg), [BubbleRaftAfter1.jpg](Gallery/BubbleRaft/BubbleRaftAfter1.jpg), [BubbleRaftAfter2.cfg](Gallery/BubbleRaft/BubbleRaftAfter2.cfg), [BubbleRaftAfter2.jpg](Gallery/BubbleRaft/BubbleRaftAfter2.jpg) as examples of input ([extended CFG format](#extended_CFG)) and output for auxiliary property coloring.

Press '[Meta](#Meta)+-','[Meta](#Meta)+=','[Meta](#Meta)+M' to change the colormaps (available: [jet](CMAP/jet.jpg), [hot](CMAP/hot.jpg), [cool](CMAP/cool.jpg), [gray](CMAP/gray.jpg), [pink](CMAP/pink.jpg), [bone](CMAP/bone.jpg), [copper](CMAP/copper.jpg), [autumn](CMAP/autumn.jpg), [spring](CMAP/spring.jpg), [winter](CMAP/winter.jpg), [summer](CMAP/summer.jpg), [hsv](CMAP/hsv.jpg)). Gradation in color only happens for atoms whose property values are between the lower threshold and the upper threshold (corresponding to colormap scales 0 and 1, respectively). Atoms whose property values are outside of the two thresholds will be rendered invisible by default. This default behavior can be toggled by pressing 'Ctrl+A', in which case those atoms whose property values are outside of the thresholds will be visible but drawn in saturated colors corresponding to colormap scales 0 and 1.

Initially, the upper threshold is taken to be the maximum of the atoms' auxiliary property values and the lower threshold is taken to be the minimum. Shift+[PageUp](#PageUp) will increase the upper threshold, Shift+[PageDown](#PageDown) will decrease the upper threshold. Ctrl+[PageUp](#PageUp) will increase the lower threshold, Ctrl+[PageDown](#PageDown) will decrease the lower threshold. This in turn controls which atoms are drawn and which atoms are invisible in the default mode. You may also use Shift+T to type in the thresholds manually. If you wish to recover the initial threshold values, press 'Ctrl+R'.

When you [reload a configuration](#reload_config) or [make a movie](#making_movie) using [animation script](#animation_script), the default threshold behavior is "floating", meaning the program reestablishes the thresholds for every new config according to the new minimum and maximum. If you do not want this behavior, pressing 'Ctrl+T' will toggle to a mode whereby the thresholds remain fixed, or "rigid".

The coloring of atoms can be restored to normal values (of chemical species) with invisibilities removed by pressing 'o'.

If you make screenshot in [eps](#making_eps) or [jpg](#making_jpg) or [png](#making_png), there will an [eps file](Gallery/BubbleRaft/BubbleRaftAfter2.jpg.cmap.eps) saved for the color scale, with numerical value labels.<a name="extra_color_patch"></a>

#### <a name="extra_color_patch">Apply Extra Color Patch:</a>

<a name="extra_color_patch"></a>

<a name="extra_color_patch">One can assign arbitrary colors and radii to atoms and bonds by putting a '.usr' patch file with the same name in the same directory as the '.cfg' file. See for example the file combination</a> [Cu.cfg](Gallery/Cu/Cu.cfg)  |  [Cu.usr](Gallery/Cu/Cu.usr). This gives ultimate control of rendering spheres and cylinders to the user. This '.usr' file should consists of lines that contain either 1, 3, 4, 5, or 6 numbers. The 6-numbers line describes bonds, the rest describes atoms.

A 6-numbers line like

<pre>0 3  0 1 0  0.2
</pre>

would mean drawing a bond between atom 0 and atom 3, with bond color "0 1 0" and radius 0.2\. Note that atom index in AtomEye starts from 0 as in C language convention (that is, the first atom in the cfg file has index 0).

A 5-numbers line like

<pre>2  0 1 0  0.7
</pre>

would mean drawing atom 2 with color "0 1 0" and radius 0.7

A 1-number line

<pre>0.7
</pre>

or a 3-numbers line

<pre>0 1 0
</pre>

or a 4-numbers line

<pre>0 1 0  0.7
</pre>

would mean drawing atom _m_ with radius 0.7 (with default color), color "0 1 0" (with default radius), or color "0 1 0" and radius 0.7, respectively. _m_ will be automatically accumulated, starting from 0.

The '.usr' patch file will be automatically loaded and applied when one presses ['](#step_config)[Insert](#Insert)'/'[Delete](#Delete)' to advance config list. See [an example of a group of files here](Gallery/CuBond/).<a name="make_atoms_invisible"></a>

#### <a name="make_atoms_invisible">Make Atoms Invisible:</a>

<a name="make_atoms_invisible"></a>

<a name="make_atoms_invisible">Ctrl+Shift+Right-click will make a certain species or coordination numbered atoms invisible. Also, when prompted to enter atom color RGB (see</a> [change atom color and/or radius](#change_atom_color_radius)), '-1 0 0' will make it invisible.

When CapsLock is at On or [Meta](#Meta) key is pressed, right-clicking on a particular atom will make it invisible.

The coloring of atoms can be restored to normal values (of chemical species) with invisibilities removed by pressing 'o'.<a name="make_bond_invisible"></a>

#### <a name="make_bond_invisible">Make Bond Invisible:</a>

<a name="make_bond_invisible"></a>

<a name="make_bond_invisible">When CapsLock is at On or</a> [Meta](#Meta) key is pressed, right-clicking on a particular bond will make it invisible.

The coloring of bonds can be restored to normal values with invisibilities removed by pressing 'o'.<a name="change_atom_color_radius"></a>

#### <a name="change_atom_color_radius">Change Atom Color and/or Radius:</a>

<a name="change_atom_color_radius">Ctrl+Shift+Left-click will prompt the user to change the RGB color and/or radius of atoms of a certain chemical species, as follows:

<pre>Change color [radius] of type-3 ("P") atoms (0.100 0.700 0.300 [1.060]):
</pre>

The first three numbers are the RGB values of</a> [Phosphorus](http://www.webelements.com/webelements/elements/text/P/key.html). The last number in the square bracket is its radius (currently taken to be the [charge radius](http://www.fhi-berlin.mpg.de/th/balsac/balm.47.html), another alternative is the [empirical atomic radius](http://www.webelements.com/webelements/properties/text/definitions/atomic-radius-emp.html)) in Angstrom. The rendered radius is the above scaled by a common factor ([scale atom radii](#scale_atom_radii)).

If you press Return after being prompted, the values would remain identical and nothing happens. If you enter one number and Return, it will be interpreted as the new radius in Angstrom. If you enter three numbers and Return, it will be interpreted as the new color. Of course, you can enter all four numbers, in which case the first three numbers will be taken as the new color.

This [link](http://alum.mit.edu/www/liju99/NetApp/JavaScript/Reader/rgb.html) may help you pick a new color. You may directly enter "255 250 250" instead of converting them to floating numbers ("1.0 0.98 0.98") first.

When CapsLock is at On or [Meta](#Meta) key is pressed, left-clicking on a particular atom will prompt the user to change its RGB color and/or radius instead of the whole species.

The coloring of atoms can be restored to normal values (of chemical species) with invisibilities removed by pressing 'o'.<a name="change_bond_color_radius"></a>

#### <a name="change_bond_color_radius">Change Bond Color and/or Radius:</a>

<a name="change_bond_color_radius"></a>

<a name="change_bond_color_radius">When CapsLock is at On or</a> [Meta](#Meta) key is pressed, left-clicking on a particular bond will prompt the user to change the RGB color and/or radius of the particular bond, as follows:

<pre>Change color [radius] of bond-1726 (0.493 0.493 0.561 [0.250]):
</pre>

The first three numbers are the RGB values of the bond. The last number in the square bracket is its radius in Angstrom.

If you press Return after being prompted, the values would remain identical and nothing happens. If you enter one number and Return, it will be interpreted as the new radius in Angstrom. If you enter three numbers and Return, it will be interpreted as the new color. Of course, you can enter all four numbers, in which case the first three numbers will be taken as the new color.

This [link](http://alum.mit.edu/www/liju99/NetApp/JavaScript/Reader/rgb.html) may help you pick a new color. You may directly enter "255 250 250" instead of converting them to floating numbers ("1.0 0.98 0.98") first.

The coloring of bonds can be restored to normal values with invisibilities removed by pressing 'o'.

As a last resort, the user can [use '.usr' file to completely redefine what bonds are drawn, their colors and radii](#extra_color_patch).<a name="shift_object_under_PBC"></a>

#### <a name="shift_object_under_PBC">Shift object under PBC:</a>

<a name="shift_object_under_PBC"></a>

<a name="shift_object_under_PBC">Shift + mouse drag to shift object under PBC. Alternatively, when CapsLock is at On, Left, Right, Up, Down, Shift+Up, Shift+Down will do equivalent things as in</a> [shift object](#shift_object), except now in PBC. Lastly, Shift+IMWheel (or Shift+Ctrl+IMWheel for quicker action) will shift object under PBC in the forward/backward direction.

Press 'z' to recover the initial PBC state, where there is no shift.<a name="change_bgcolor"></a>

#### <a name="change_bgcolor">Change Background Color:</a>

<a name="change_bgcolor">

Press 'd', xterm will pop out with inquiry,

<pre>Change background color (0.000 0.000 0.000):
</pre>

and you should input three real numbers from 0 to 1 for the RGB values. If you press Return, the numbers in the brackets will be taken as the input. They are just the current background colors, so nothing will change.</a><a name="change_total_view_angle">

#### Change View Angle Amplification:

</a>

<a name="change_total_view_angle">Shift+Home/End changes the total view angle that the viewport spans. The smaller it is, the larger the image appears on screen. When the session starts, the total span of view is 60 degrees.</a><a name="new_clone_quit_viewport"></a>

#### <a name="new_clone_quit_viewport">New, Clone, Quit Viewport:</a>

<a name="new_clone_quit_viewport"></a>

<a name="new_clone_quit_viewport">Occasionally you want another perspective on the same configuration. Press 'F4' to new a viewport where all parameters take default values. Press 'c' to clone a viewport that is exactly the same as the current one. Press 'q' to quit a certain viewport and free its memory.</a><a name="system_status"></a>

#### <a name="system_status">Print System Status:</a>

<a name="system_status"></a>

<a name="system_status">Press 's' to print out system status.</a><a name="find_atom"></a>

#### <a name="find_atom">Find an Atom:</a>

<a name="find_atom"></a>

<a name="find_atom">Press 'f' to locate an atom by entering its index (0 to number_of_particles-1). The</a> [anchor](#anchor_control) will then be set on that atom, and by [pulling closer/away](#pull_closer_away) it will be quickly obvious which atom it is.<a name="goto_position"></a>

#### <a name="goto_position">Go to Position:</a>

<a name="goto_position"></a>

<a name="goto_position">Press 'g' to go to certain position, maybe somewhere inside the box, by entering three reduced coordinates.</a><a name="resize_window"></a>

#### <a name="resize_window">Resize Window</a>

<a name="resize_window"></a>

<a name="resize_window">Just dragging the corner of the window will do, but for more precise control using terminal input, press 'Ctrl+S'. The recommended sizes are 320x240, 640x480, 800x600, 1024x1024 ... for certain video compressors.</a><a name="making_jpg"></a>

#### <a name="making_jpg">Making</a> [jpeg](http://www.jpeg.org/) Screenshot:

Press 'j' to make .jpg screenshot.<a name="making_png"></a>

#### <a name="making_png">Making</a> [png](http://libpng.org/pub/png/) Screenshot:

Press 'p' to make .png screenshot.<a name="making_eps"></a>

#### <a name="making_eps">Making Encapsulated PostScript Screenshot:</a>

<a name="making_eps">

Press 'e' to make high-resolution .eps screenshot (try scale & print).

</a>

<a name="making_eps"></a>
*   <a name="making_eps"></a>[DNA.eps](Gallery/DNA/DNA.eps)
*   [SiC_NanoXtal.eps](Gallery/SiC_NanoXtal/SiC_NanoXtal.eps)

<a name="toggle_viewer">

#### Toggle Auto-Invoking Shell Viewer:

</a>

<a name="toggle_viewer">Press 'v' to toggle auto-invoking shell viewer for screenshots. '</a>[xv](http://www.trilon.com/xv/)', '[ghostview](http://www.cs.wisc.edu/~ghost/)', '[gv](http://wwwthep.physik.uni-mainz.de/~plass/gv/)' are some of the applications AtomEye tries to find under the current PATH.<a name="load_new_config"></a>

#### <a name="load_new_config">Load New Config:</a>

<a name="load_new_config"></a>

<a name="load_new_config">Press 'F9' to load in new config file.</a><a name="reload_config"></a>

#### <a name="reload_config">Reload Configuration:</a>

<a name="reload_config"></a>

<a name="reload_config">Press 'F10' to reload the config file, in case it has been refreshed in the meantime.</a><a name="step_config"></a>

#### <a name="step_config">Sequential Config List Browsing:</a>

<a name="step_config">

AtomEye tries to compile a list of config files sequentially similar to the current config. For example, if we are viewing "a005.cfg", AtomEye will automatically line up "a000.cfg", "a001.cfg", .. in front "a005.cfg", and "a006.cfg", "a007.cfg", ... after "a005.cfg", if these files indeed exist.

</a>

<a name="step_config">Press '</a>[Insert](#Insert)' to backtrack config list and '[Delete](#Delete)' to advance config list. Press 'Ctrl+[Insert](#Insert)' to go to the first of the config list and 'Ctrl+[Delete](#Delete)' to go to the last of the config list. Loop-back is supported at the two terminations. Press 'Ctrl+F12' to change the stepping of advance/backtrack; the default stepping is 1.<a name="color_tiling"></a>

#### <a name="color_tiling">Define and Trace Color Tiling Blocks:</a>

<a name="color_tiling"></a>

<a name="color_tiling">Press 'F2' to (re)define color tiling blocks. Then, sometime later on (after a new config has been loaded, see for instance</a> [Sequential Config List Browsing](#step_config)), press 'F3' to show the atoms in previously allocated block colors - be mindful, however, that the later configurations must have the same number of atoms as the previous configuration when the color blocks are defined. If you no longer need to trace anymore, you can free this extra piece of coloring memory by 'Ctrl+F2'.<a name="cutting_planes"></a>

#### <a name="cutting_planes">Creating and Manipulating Cutting Planes</a>

<a name="cutting_planes"></a>

<a name="cutting_planes">'Shift+[0-9,a-f]' will toggle one of the 16 available cutting planes. The most recently activated cutting plane also gains the focus. The focused cutting plane can be advanced/retracted by 'Shift+RightArrow/LeftArrow', its sense can be flipped by 'Shift+P', and it can be deleted from memory by 'Shift+</a>[Meta](#Meta)+[0-9,a-f]'.

A cutting planes is created when being toggled for the first time. The user will be asked to supply six numbers: dx,dy,dz,s0,s1,s2, where (dx,dy,dz) is the normal vector (doesn't have to normalized) of the cutting plane in Cartesian frame, (s0,s1,s2) is a point on the plane in reduced coordinates from 0 to 1\. By default dx=1,dy=1,dz=1, which means the 111 plane, and (s0,s1,s2) is the current [anchor](#anchor_control) position. Usually the user just needs to input three numbers dx,dy,dz. There are three ways to input dx,dy,dz,s0,s1,s2:

*   Manually compute them, "Rain Man" style.
*   Right-click on atom A once, then right-click on atom B _twice_, then 'Shift+[0-9,a-f]' and press return (accept the defaults). This will create a cutting plane between atom A and B, keeping atom B.
*   Right-click on atom A once, atom B once, atom C once, then 'Shift+[0-9,a-f]' and press return (accept the defaults). This will create a cutting plane containing A-B-C, with right-handed normal.

You can shift an activated cutting plane to the current anchor position by 'Shift+Ctrl+[0-9,a-f]'.

Besides filtering the atoms, the cutting planes are also represented by wireframes of their sections with the H-box. A focused cutting plane's wireframe's color can toggled down to invisible by 'Shift+I'.<a name="save_atom_indices"></a>

#### <a name="save_atom_indices">Save Atom Indices in a File</a>

<a name="save_atom_indices"></a>

<a name="save_atom_indices">To create a dislocation or a crack at the desired place and inclination, it is usually the easiest to identity two adjacent crystallographic planes, and apply a force or displacement dipole. To identify the atoms involved in this operation, first create a</a> [cutting plane](#cutting_planes) by for instance pressing 'Shift+0' on [FCC10x10x10.cfg](mul/FCC10x10x10.cfg). Then, you need to click on three atoms sequentially on the exposed surface (to visualize the atoms you click, you may have CapsLock 'on'), and then press ';'

You will be prompted:

<pre>Interplanar spacing [A]  z-margin [A]  xy-tolerance (2.5 0.01 0.01):
</pre>

The defaults are generally OK for selecting two atomic planes in fcc Cu with (111) interplanar spacing 2.0871 Angstrom. Press return and you will then see:

<pre>down=4 and up=1 atoms selected in filter.
Save the selected atoms to (FCC10x10x10.idx):
selected atom indices [0-31999] saved to FCC10x10x10.idx.
</pre>

which means 4 atoms in the "down" plane and 1 atom in the "up" plane are selected, and the atom indices (remember that AtomEye always counts from 0) are saved in "FCC10x10x10.idx", which looks like:

<pre>% more FCC10x10x10.idx
4 1
14267
14550
14576
17421
17760
</pre>

The first 4 indices are atoms on the "down" parallelogram, followed by the 1 atom on the "up" parallelogram. Note that the "up" direction is the right-handed rotational axis represented by the three sequentially clicked atomic positions. In other words, you must click on the "down" atoms in anti-clockwise fashion. The make sure, you can explicitly check the printout

<pre>"up" is [0.57735 0.57735 0.57735]
check it agrees with your mirror normal...
</pre>

with your desired "up" direction.

Your favorite MD or conjugate gradient code can then load in this "FCC10x10x10.idx" file and have some fun moving these selected atoms. Before you quit, you may identify the Burgers vector as well by [inquiring geometrical info](#inquire_geometrical_info).<a name="making_movie"></a>

#### <a name="making_movie">Making movie:</a>

<a name="making_movie"></a>

<a name="making_movie"></a>
2.  <a name="making_movie">One must first have a sequence of configuration files, say "00001.cfg", "00002.cfg", ..., "00005.cfg". (a sample program</a> [sequence.f](sequence.f) shows how to do it in FORTRAN).<a name="animation_script"></a>
3.  <a name="animation_script">AtomEye can render a series of .</a>[jpg](#making_jpg), .[png](#making_png), or .[eps](#making_eps) images according to a script (default name = "scr_anim"), which looks like:

<pre>90
00001.cfg Pic/00001.jpg
00002.cfg Pic/00002.jpg
00003.cfg Pic/00003.jpg
00004.cfg Pic/00004.jpg
00005.cfg Pic/00005.jpg
</pre>

The first integer (from 0 to 100) specifies the quality of the images. The larger it is, the less the compression and the better the quality. Usually 90 is very good and 80 is good enough. Rows of string pairs then follow. The first string is the input configuration filename. The second string is the output image filename. Depending on the suffix, .[jpg](#making_jpg), .[png](#making_png), or .[eps](#making_eps) screenshots will be saved. One can also use absolute pathnames. It is better to separate the configuration files with the image files, therefore "Pic/00001.jpg" and so on is recommended. "[scr_anim](scr_anim)" itself can be created for example by using the following Matlab script [create_scr_anim.m](create_scr_anim.m). If "[scr_anim](scr_anim)" does not exist, AtomEye will offer to create a default one, in which case .[jpg](#making_jpg) format is used with quality 90, and the frames will go to the [Jpgs/](Compression/Jpgs/) directory (created by AtomEye if it doesn't exist already).
4.  One now needs to pick a view point and other rendering options for the movie. To do that, run "A 00005.cfg &", [choose screen size](#resize_window) (screen sizes like 256 by 256, 512 by 512, 768 by 512 are preferred by the video compressors), angle and coloring etc., then press 'y'. AtomEye asks you which animation script to use (default = "scr_anim"), and will start to render "00001.cfg", "00002.cfg", ..., "00005.cfg" sequence using the chosen viewpoint and other options.
5.  On Linux, you can use the powerful [MPlayer / MEncoder](http://www.mplayerhq.hu/) utility to compress sequential images to a movie. See an example shell script [vidcompress](Compression/vidcompress) and its [result](Compression/). MPEG4 codec is recommend for portability. You may also use commercial software like [Animation Shop](http://www.jasc.com/) or [Adobe Premiere](http://www.adobe.com/) to do the same thing. Use [VirtualDub](http://www.virtualdub.org/) to do extra editing.
6.  The AVI or MPEG4 movie should be directly playable on Microsoft Windows Media Player. It may also be directly inserted into Microsoft PowerPoint (pull-down menu: Insert / Movies and Sounds / Movie from File). [MPlayer](http://www.mplayerhq.hu/) is an excellent free AVI and MPEG4 movie player on Linux.

<a name="Meta">

#### What is the Meta key?

Some environments forbid the use of Meta or Alt key. In those cases, having CapsLock at "on" works as if the Meta key is pressed.

* * *

</a><a name="MacKeybinding">

## Mac Keybinding

</a><a name="PageUp">PageUp (Windows) = fn + uparrow (Mac)</a>

<a name="PageUp"></a><a name="PageDown">PageDown (Windows) = fn + downarrow (Mac)</a>

<a name="Delete">Delete (Windows) = fn + delete (Mac)</a>

<a name="Insert">Insert (Windows) = shift + fn + delete (Mac), or</a>

<a name="Insert">Insert (Windows) = \ (Mac)</a><a name="download"></a>

## <a name="download">Binary Release</a>

<a name="download">These are raw binaries. Right-click on the link and "Save Target As..." to one of your directories. Then run "chmod 755" on the file. To</a> [test](#usage), you need to save a [CFG](#standard_CFG) file as well, such as [cnt8x3.cfg](Gallery/Nanotube/cnt8x3.cfg).

*   [i686 Linux](Download/A.i686)
*   [Alpha Linux GLIBC2.1](Download/A.alpha)
*   [Sgi Irix](Download/A.IRIX)
*   [Sgi Irix64](Download/A.IRIX64)
*   [Sun Solaris](Download/A.Sun)
*   [HP UX](Download/A.HPUX)
*   [Windows](Download/A3-20130625.exe) with [Cygwin/X](http://x.cygwin.com/) ([README.txt](Download/CYGWIN/README.txt))
*   [Alpha Tru64 UNIX](Download/A.OSF1)
*   Mac OS X ([PowerPC](Download/A.Darwin), [Intel386](Download/A.DarwinIntel386), [Intel x86_64](Download/A.DarwinIntelx86_64)) with [Darwin](http://developer.apple.com/darwin/) (see [A](http://www.xfree86.org/4.3.0/Darwin2.html), [B](http://www.mrcla.com/XonX/FAQ.html#3buttons), [C](http://www.xdarwin.org/faq/#fakebuttons) for button issues)

#### When I get the chance, I try to make up-to-date executables for different platforms, but only the i686 Linux one is <u>guaranteed</u> to be of the latest stable version with full functionality.

* * *

<a name="FAQ">

## Frequently Asked Questions

</a>

<a name="FAQ">

1. Some Linux machines give semi-transparent Atomeye windows, how to correct it?

**Answer**: Try adding

<pre>export XLIB_SKIP_ARGB_VISUALS=1
</pre>

to your .bashrc file.

2. Is there a way to display what value each color in a color-mappped variable represent (color bar)?

**Answer**: If you save the screenshot by pressing 'j' or 'p' or 'e', there should be an extra colorbar file saved in .eps with numerical labels.


3. I encounter the following error message:

<pre>ATOM_COORDINATION_MAX = 24 exceeded</pre>

or
<pre>error: Imakespace: min=0 max=2744 jammed between i=1262 and i+1.</pre>

and AtomEye quits. What's the problem?

**Answer**: The above error message means some atoms in the configuration are getting too close to each other. The number of atoms within the default cutoff radii of first-nearest-neighbors exceeds 24\. This usually means there is some pathology in the configuration (maybe you have miscalculated the atomic geometry? maybe your time-integrator has blown up?).

</a>

<a name="FAQ">To see what is in the configuration, artificially scale up your supercell 10 times by adding</a> [one optional line](#standard_CFG_example) after the 'Number of particles = ' line:

<pre>A = 10 Angstrom (basic length-scale)
</pre>

Most likely this would reduce all atoms to coordination-0\. But at least you can now load your configuration into AtomEye, see the atoms, and debug the configuration.

* * *

<a name="bugreport">

## Bug report

There will always be bugs in AtomEye. When you think you have encountered one, please do the following:</a>

<a name="bugreport">*   To the best of your capability, try to make the bug repeatable.</a>
*   <a name="bugreport">Write</a> [me](http://164.107.79.177/lij/card.html) an email, with detailed, step-by-step instructions, so I can reproduce the bug; and attach all necessary files in the email.

* * *

<a name="history">

## History

I started using</a> [Rasmol](http://www.bernstein-plus-sons.com/software/rasmol/) in 1997 and was quite impressed. But its as well as [PDB](http://www.rcsb.org/pdb/) file format's limitations for large-scale molecular dynamics (MD) simulations became painfully obvious after some time. My colleagues Dongyi Liao, Wei Cai and myself all started developing our own molecular visualization codes. It was a very friendly and helpful environment even though we differ in coding styles, and a lot of ideas were exchanged. Some of my earlier attempts are placed [here](../).

Gradually it was realized that [X-Window](http://www.x.org/) does not provide sufficient native capability. There are two options. One is to program in [OpenGL](http://www.opengl.org/), the other is to basically rewrite the graphics functionality of Xlib starting from scratch (antialised lines, clipping, tilted ellipses, etc.) and just to use [X-Window](http://www.x.org/) as a pipe (MIT-SHM extension / XPutImage). I took the second approach based on the observation that one needs to render _very_ few types of objects in massive quantities: spheres (as atoms), cylinders (as bonds), and points (as charge density). Therefore very efficient routines, such as sphere caching, can be written, which is included in the _libAX_ (Advanced X) library. [Graphics card](http://www.nvidia.com) has evolved, but so has CPU / main memory, as they essentially are from the same technology. Therefore in foreseeable future, it is highly unlikely that this approach loses to [OpenGL](http://www.opengl.org/) which depends on polygon assemblage.

In addition to _libAX_, a crucial element of AtomEye is the order-N treatment of atomistic configurations for bond calculation and coordination number / local strain coloring of atoms. This library, _libAtoms_, is shared between AtomEye and an order-N MD code. We assume always that the configuration is under periodic boundary condition (PBC) in a parallelepiped box. The rationale is that it is not too difficult to represent a non-PBC configuration as a PBC configuration, but not the other way around.

AtomEye was first compiled in 1999, the executable is called [A](#download). The early users are myself, Dongyi Liao, Wei Cai, Dongsheng Xu, Jinpeng Chang and Shigenobu Ogata. Extensive debugging and improvements have been made to the code, driven by user requests. AtomEye's screenshots have appeared on the cover of [_Nature Materials_](Gallery/Covers/NatureMaterialsNov2007Volume6No11/IE.pdf), [_PNAS_](Gallery/Covers/PNASFeb2007Volume104No9/InterfacialPlasticity.pdf) and [_Physical Review Letters_](Gallery/Covers/PRLJan2008Volume100No2/cv100002.pdf).

* * *

<address></address>

![](http://li.mit.edu/cgi-bin/AtomEyeCounter.cgi)
