import os
import json
import subprocess
from pathlib import Path

stl_log_filename = '.stl.log'

def parse_json(json_file):
    with open(json_file, 'r') as f:
        s = f.read()
        obj = json.loads(s)

    l = obj['children'][0]['params']['length']
    w = obj['children'][0]['params']['width']
    h = obj['children'][0]['params']['height']
    n = obj['children'][0]['children'][0]['children'][0]['params']['noseLength']
    r = obj['children'][0]['children'][1]['children'][0]['children'][0]['params']['radius']
    t = obj['children'][0]['children'][2]['children'][0]['params']['tailLength']
    e = obj['children'][0]['children'][2]['children'][0]['params']['endRadius']
    
    return [l, w, h, n, r, t, e]

def gen_stl(length, width, height, noseLength, radius, tailLength, endRadius, name=None):
    # Generate boxfish stl file
    data = {
            "type": "root",
            "params": {"segments":128},
            "children": [
                {
                    "type": "SingleHull",
                    "params": {
                        "length": length, #2400,
                        "width": width, #700,
                        "height": height, #700
                        },
                    "children": [
                        {
                            "type": "Front",
                            "children": [
                                {
                                    "type": "LinearNoseCone",
                                    "params": {
                                        "noseLength": noseLength, #15
                                        }
                                    }
                                ]
                            },
                        {
                            "type": "Midsection",
                            "children": [
                                {
                                    "type": "Midshape",
                                    "children": [
                                        {
                                            "type": "Boxfish",
                                            "params": 
                                            {
                                                "rounded": 1,
                                                "radius": radius, #200
                                                }
                                            }
                                        ]
                                    }
                                ]
                            },
                        {
                            "type": "Tail",
                            "children": [
                                {
                                    "type": "LinearTailCone",
                                    "params": {
                                        "tailLength": tailLength, #15,
                                        "endRadius": endRadius, #100
                                        }
                                    }
                                ]
                            }
                        ]
                    }
                ]
            }

    # stl_filename = f'{stl_dir}/boxfish_{length}_{width}_{height}_{noseLength}_{radius}_{tailLength}_{endRadius}.stl'
    if name is None: name = 'tmp'
    stl_filename = f'{name}.stl'
    json_str = json.dumps(data, separators=(',',':'))
    with open(f'{name}.json', 'w') as f:
        f.write(json_str)
    stl_cmd = f"npx jscad jscad-generator.js -o {stl_filename} --jsonComponentModel {json_str}"
    with open(stl_log_filename, 'w') as outfile:
        subprocess.run(stl_cmd.split(' '), stdout=outfile, stderr=outfile)
    return stl_filename

def parse_volume():
    with open(stl_log_filename, 'r') as f:
        lines = f.readlines()
        for line in list(reversed(lines)):
            if 'area' in line and 'volume' in line:
                return json.loads(line)['volume']

def parse_area():
    with open(stl_log_filename, 'r') as f:
        lines = f.readlines()
        for line in list(reversed(lines)):
            if 'area' in line and 'volume' in line:
                return json.loads(line)['area']

def gen_designs(sol, run_name=""):
    path = Path(f"finals/{run_name}")
    path.mkdir(parents=True, exist_ok=True)
    fairing_name = "fairing"
    payload_name = "OBS_V1_Disp"
    pv_name = "Float_Brick_1"

    gen_stl(*(sol.params), str(path / fairing_name))

    v_length = sol.params[0] / 1e3
    v_width = sol.params[1] / 1e3
    v_height = sol.params[2] / 1e3

    obj = {}
    for item in sol.pack.packer.bins[0].items:
        name = item.name
        position = [float(pos) / 1e3 for pos in item.position]
        rotation = item.rotation_type
        dimension = [float(dim) / 1e3 for dim in item.get_dimension()]
        obj[name] = {
                'position': position,
                'rotation': rotation,
                'dimension': dimension,
                }

    with open(str(path / f"{pv_name}.csv"), 'w') as f:
        f.write('"Float Brick,  Definition Sheet",,,\n')
        f.write("Material - Syntactic Foam (4000 meter rating),,,\n")
        f.write("brick_width (m),brick_length (m),brick_height  (m),Density of Material (kg/m^3)\n")
        l, w, h = obj['pv']['dimension']
        f.write(f"{l},{w},{h},{0}\n")

    with open(str(path / 'Vehicle_Def.csv'), 'w') as f:
        f.write("Vehicle Definition File,,,,,,,,,\n")
        f.write("Water Density (kg/m^3),1027,,,,Module Position,,,,\n")
        f.write('Module File (.STEP),Module Density (kg/m^3),Module_Disp file (.STEP),Module_Disp Density (kg/m^3),"Module Type (WET, DRY or PV)",X (m),Y (m),Z (m),"Axis of Rotation (X, Y, Z)",Rotation Angle (Deg.)\n')

        f.write(f"{fairing_name},0,{fairing_name},0,WET,0,0,0,X,0\n")
        x, y, z = obj['payload']['position']
        f.write(f"{payload_name},0,{payload_name},0,WET,{x},{y},{z},X,0\n")
        x, y, z = obj['pv']['position']
        f.write(f"{pv_name},0,{pv_name},0,WET,{x},{y},{z},X,0\n")

    files = ['Assembly_macro_1.FCMacro', 'OBS_V1_Disp.STEP', 'Float_Brick.FCMacro', 'Start_Macro_3.FCMacro']
    # for filename in files:
    #     cmd = f"ln -s ./cad/{filename} ./{str(path / filename)}"
    #     os.system(cmd)
    cmd = f"cp -r ./cad/ ./{str(path)}"
    os.system(cmd)

    print(f"> Design generated. See {str(path)} for related files.")
    
