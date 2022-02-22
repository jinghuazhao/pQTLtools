#!/usr/bin/env python
from __future__ import print_function
import sys
import os

args = sys.argv
print('')

help_message = """
--------------------------
STRING node color changer.

USAGE:

    (1) Download the SVG_FILE of the network you want to modify.
        You can find it in the EXPORTS tab below the picture of the network. 

    (2) Generate node COLOR_TABLE based on this SVG file:

        $ python change_STRING_colors.py -s SVG_FILE
  
        This command will output COLOR_TABLE file (color_table.tsv) which you can find in
        the same directory as the script.

    (3) Modify the colors in the COLOR_TABLE (color_table.tsv) in Excel, Google Docs or any simple text editor.
        Save the modify sheet as a TSV or plain text file. 
        
        As colors you can use use all the standard HTML color names:
           - red (see: https://htmlcolorcodes.com/color-names/)
        ...or RGB values in this syntax:
           - rgb(255,0,0)
        ...or hex color values:
           - #FF0000
        ...or if you remove the node completely from the file, it will turn white.


    (4) Generate new SVG file with the modified colors:
 
        $ python change_STRING_colors.py -s SVG_FILE -c COLOR_TABLE

        The command will generate new SVG file in the same folder as the original SVG file.
        To generate the raster graphics open the modified SVG in any program that
        can read vector graphics (Adobe Illustrator, Affinity Designer or Inkscape (free) and
        export it to PNG, BMP or JPG format.  



""" 

insert_mode = True

if "-s" in args:

        i = args.index("-s")
        
        try:
            svg_file = args[i+1]
        except Exception as e:
             print("ERROR: SVG file not specified")
             sys.exit(help_message)

        if not os.path.exists(svg_file):
            print("ERROR: SVG file not found")
            sys.exit(help_message)
else:
   print("ERROR: -s option not specified")
   sys.exit(help_message)

if "-c" in args:
        i = args.index("-c")
        try:
            color_file = args[i+1]
        except Exception as e:
             print("ERROR: color file not specified")
             sys.exit(help_message)

else:
   print("Generating color table for the SVG (color_table.tsv)")
   insert_mode = False
   color_file = "color_table.tsv"

nodes_color = {}

if insert_mode:
    
    if not os.path.exists(color_file):
        sys.exit("ERROR: color table file not found...")

    print("Reading colors...")

    try:
        for line in open(color_file):
    
            l = line.strip().split("\t")
            if len(l) >= 2:
                orig_node = l[0]
                node = l[0].lower()
                color = l[1].lower()
                color = color.replace(" ", "")
                color = color.replace("_", "")
    
                if node not in nodes_color:
                    nodes_color[node] = color
                else:
                    print("WARNING: Node duplicated in the color file  %s" % orig_node)
    
    except Exception as e:
        print("ERROR parsing color file")
        sys.exit(help_message)
    
    print("Numbers of nodes read from the file:", len(nodes_color))

if not insert_mode:
    fh_out = open(color_file, "w")
else:
    new_svg_file = svg_file.rsplit(".", 1)[0] + ".new_colors.svg"
    fh_out = open(new_svg_file, "w")

if insert_mode:
    print("Parsing SVG file...")
else:
    print("Parsing SVG file and writing colors_table.tsv...")

assigned_colors = []
node = ""
for line in open(svg_file):

    line = line.replace("'", '"')

    if 'class="nwbubblecoloredcircle' in line:
        color = line.split("fill=\"")[-1].split()[0].strip('"')

    if '<text fill="black"' in line or '<text fill="rgb(0,0,0)"' in line:
        node = line.split('">')[-1].split("</")[0]
        orig_node = node
        node = node.lower()
        if not insert_mode:
            print(orig_node, color, sep="\t", file=fh_out)
            assigned_colors.append((node, color))

        elif node in nodes_color:
            assigned_colors.append((node, nodes_color[node]))
        else:
            assigned_colors.append((node, 'rgb(255,255,255)'))

if not insert_mode:
    if len(assigned_colors):
        print("Number of nodes found in the SVG file: %s" % len(assigned_colors))
    else:
        print("ERROR: No nodes found in the SVG ")
        exit(0) 

if insert_mode:

    for line in open(svg_file):
    
        line = line.replace("'", '"')

        if 'class="nwbubblecoloredcircle' in line:
            color = line.split('fill="')[-1].split()[0].strip('"')
            l = line.split('fill="')
            lcolor = l[-1].split()
            node, assigned_color = assigned_colors.pop(0)
            l[-1] = 'fill="%s" %s\n' % (assigned_color, " ".join(lcolor[1:]))
            line = " ".join(l)
    
        fh_out.write(line)

fh_out.close()

if insert_mode:
    print("Recolored SVG file save as: " + new_svg_file)
else:
    print("Color table saved as: color_table.tsv")

print("All done.")
