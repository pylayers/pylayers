
# PyLayers Overview

## Geographical Information 

+ [Description of Building Layouts](1-GIS/Layout.html)
+ [Editing a Layout](1-GIS/LayoutEditor.html)
+ [Effects of changing constitutive properties of a
Layout](1-GIS/Multisubsegments.html)
+ [Handling of Geographical Data: The Earth Zone class](1-GIS/Ezone.html)

## Antenna and Propagation

### Antennas

+ [Description of Antennas](2-AP/Antenna.html)
+ [Scalar Spherical Harmonics Representation](2-AP/AntennaSSH.html)
+ [Vector Spherical Harmonics Representation](2-AP/AntennaVSH.html)

### Slabs and Materials


+ [SlabsMaterials](2-AP/SlabsMaterials.html)


### MultiWall model

+ [Coverage](2-AP/Coverage.html)
+ [MultiWall](2-AP/MultiWall.html)

### Signatures and Rays

+ [The Ray Signatures](2-AP/Signatures.html)

+ [The transmission Channel](2-AP/Channel.html)

## Handling time domain signals

+ [The Bsignal Class](3-PHY/Bsignal.html)

## Human Mobility

+ [Simulnet Configuation](5-SIM/SimulNetConfig.html)
+ [Large scale mobility `simulnet`](4-MOB/Mobility.html)

+ [Body scale mobility ](4-MOB/Body.html)

## Ray Tracing Simulation Examples

+ [Link Simulations](5-SIM/LinkSimulation.html)
+ [DLR WHERE2 Scenario](5-SIM/DLR-WHERE2.html)
+ [PTIN WHERE2 Scenario](5-SIM/PTIN.html)
+ [WHERE1 - M1 Scenario](5-SIM/WHERE1-M1.html)
+ [WHERE1 - M1 Scenario (bis)](5-SIM/Where1M1.html)
+ [Aggregated CIR](5-SIM/AggregatedCIR.html)

## PyLayers Classes and Tools

+ [The Cone class](8-MISC/Cone.html)
+ [Geometrical Utility Functions](8-MISC/Geomutil.html)


    from IPython.core.display import HTML
    
    def css_styling():
        styles = open("../styles/custom.css", "r").read()
        return HTML(styles)
    css_styling()




<style>
    @font-face {
        font-family: "Computer Modern";
        src: url('http://mirrors.ctan.org/fonts/cm-unicode/fonts/otf/cmunss.otf');
    }
    div.cell{
        width:800px;
        margin-left:16% !important;
        margin-right:auto;
    }
    h1 {
        font-family: Helvetica, serif;
    }
    h4{
        margin-top:12px;
        margin-bottom: 3px;
       }
    div.text_cell_render{
        font-family: Computer Modern, "Helvetica Neue", Arial, Helvetica, Geneva, sans-serif;
        line-height: 145%;
        font-size: 130%;
        width:800px;
        margin-left:auto;
        margin-right:auto;
    }
    .CodeMirror{
            font-family: "Source Code Pro", source-code-pro,Consolas, monospace;
    }
    .prompt{
        display: None;
    }
    .text_cell_render h5 {
        font-weight: 300;
        font-size: 22pt;
        color: #4057A1;
        font-style: italic;
        margin-bottom: .5em;
        margin-top: 0.5em;
        display: block;
    }
    
    .warning{
        color: rgb( 240, 20, 20 )
        }  
</style>
<script>
    MathJax.Hub.Config({
                        TeX: {
                           extensions: ["AMSmath.js"]
                           },
                tex2jax: {
                    inlineMath: [ ['$','$'], ["\\(","\\)"] ],
                    displayMath: [ ['$$','$$'], ["\\[","\\]"] ]
                },
                displayAlign: 'center', // Change this to 'center' to center equations.
                "HTML-CSS": {
                    styles: {'.MathJax_Display': {"margin": 4}}
                }
        });
</script>



