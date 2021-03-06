# GUI for DMS project, powered by TKinter
# Created by Kai, 2018
import matplotlib
matplotlib.use('TkAgg')
import tkinter as tk
from tkinter import ttk
from tkinter.filedialog import askopenfilename
from pathlib import Path
from ast import literal_eval as make_tuple
import wrapper
#from wrapper import *

class App:
    def __init__(self, master):
        master.title("GUI wrapper")
        master.resizable(False, False)
        master.geometry("520x510") # If using this, must import matplotlib
        self.master = master
        self.Dict = dict({'wt_file_button': 'placeholder'})

        self.wt_file_label  = ttk.Label(self.master, text='Select WT seq file: ')
        self.wt_file_label.grid(row=0, column=0)
        self.wt_file_button = ttk.Button(self.master, text='Click to choose')
        self.wt_file_button.config(command=lambda: self.openfile(self.wt_file_button, "wt_file_button"))
        self.wt_file_button.grid(row=0, column=2)

        self.exp_0_label = ttk.Label(self.master, text='Experiemnt 0: \n(input)')
        self.exp_0_label.grid(row=1, column=0, rowspan=2)
        self.exp_0_name  = ttk.Label(self.master, text='Assay name: ')
        self.exp_0_name.grid(row=1, column=1)
        self.exp_0_name_entry  = ttk.Entry(self.master, width=14)
        self.exp_0_name_entry.grid(row=1, column=2)
        self.exp_0_name_button = ttk.Button(self.master, text='confirm')
        self.exp_0_name_button.grid(row=1, column=3)
        self.exp_0_name_button.config(command=lambda: self.on_button(self.exp_0_name_button,
                                                                     self.exp_0_name_entry, 'exp_0_name_entry'))

        self.exp_0_name_button.grid(row=1, column=3)
        self.exp_0_file1_label = ttk.Label(self.master, text='Sequence file1: ')
        self.exp_0_file1_label.grid(row=2, column=1)
        self.exp_0_file1_button = ttk.Button(self.master, text='Click to choose')
        self.exp_0_file1_button.grid(row=2, column=2)
        self.exp_0_file1_button.config(command=lambda: self.openfile(self.exp_0_file1_button, "exp_0_file1_button"))
        self.exp_0_file2_label = ttk.Label(self.master, text='Sequence file2: ')
        self.exp_0_file2_label.grid(row=3, column=1)
        self.exp_0_file2_button = ttk.Button(self.master, text='Click to choose')
        self.exp_0_file2_button.grid(row=3, column=2)
        self.exp_0_file2_button.config(command=lambda: self.openfile(self.exp_0_file2_button, "exp_0_file2_button"))

        self.exp_1_label = ttk.Label(self.master, text='Experiemnt 1: ')
        self.exp_1_label.grid(row=4, column=0)
        self.exp_1_name  = ttk.Label(self.master, text='Assay name: ')
        self.exp_1_name.grid(row=4, column=1)
        self.exp_1_name_entry  = ttk.Entry(self.master, width=14)
        self.exp_1_name_entry.grid(row=4, column=2)
        self.exp_1_name_button = ttk.Button(self.master, text='confirm')
        self.exp_1_name_button.grid(row=4, column=3)
        self.exp_1_name_button.config(command=lambda: self.on_button(self.exp_1_name_button,
                                                                     self.exp_1_name_entry, 'exp_1_name_entry'))

        self.exp_1_file1_label = ttk.Label(self.master, text='Sequence file1: ')
        self.exp_1_file1_label.grid(row=5, column=1)
        self.exp_1_file1_button = ttk.Button(self.master, text='Click to choose')
        self.exp_1_file1_button.grid(row=5, column=2)
        self.exp_1_file1_button.config(command=lambda: self.openfile(self.exp_1_file1_button, "exp_1_file1_button"))
        self.exp_1_file2_label = ttk.Label(self.master, text='Sequence file2: ')
        self.exp_1_file2_label.grid(row=6, column=1)
        self.exp_1_file2_button = ttk.Button(self.master, text='Click to choose')
        self.exp_1_file2_button.grid(row=6, column=2)
        self.exp_1_file2_button.config(command=lambda: self.openfile(self.exp_1_file2_button, "exp_1_file2_button"))

        self.exp_2_label = ttk.Label(self.master, text='Experiemnt 2: ')
        self.exp_2_label.grid(row=7, column=0)
        self.exp_2_name  = ttk.Label(self.master, text='Assay name: ')
        self.exp_2_name.grid(row=7, column=1)
        self.exp_2_name_entry  = ttk.Entry(self.master, width=14)
        self.exp_2_name_entry.grid(row=7, column=2)
        self.exp_2_name_button = ttk.Button(self.master, text='confirm')
        self.exp_2_name_button.grid(row=7, column=3)
        self.exp_2_name_button.config(command=lambda: self.on_button(self.exp_2_name_button,
                                                                     self.exp_2_name_entry, 'exp_2_name_entry'))

        self.exp_2_file1_label = ttk.Label(self.master, text='Sequence file1: ')
        self.exp_2_file1_label.grid(row=8, column=1)
        self.exp_2_file1_button = ttk.Button(self.master, text='Click to choose')
        self.exp_2_file1_button.grid(row=8, column=2)
        self.exp_2_file1_button.config(command=lambda: self.openfile(self.exp_2_file1_button, "exp_2_file1_button"))
        self.exp_2_file2_label = ttk.Label(self.master, text='Sequence file2: ')
        self.exp_2_file2_label.grid(row=9, column=1)
        self.exp_2_file2_button = ttk.Button(self.master, text='Click to choose')
        self.exp_2_file2_button.grid(row=9, column=2)
        self.exp_2_file2_button.config(command=lambda: self.openfile(self.exp_2_file2_button, "exp_2_file2_button"))

        self.amplicon_1_label = ttk.Label(self.master, text='Amplicon1 range:\nformat: 1-XXX', width=15)
        self.amplicon_1_label.grid(row=10, column=0, rowspan=2)
        self.amplicon_1_entry = ttk.Entry(self.master, width=12)
        self.amplicon_1_entry.grid(row=10, column=1, rowspan=2)
        self.amplicon_1_button = ttk.Button(self.master, text = 'confirm')
        self.amplicon_1_button.grid(row=10, column=2, rowspan=2)
        self.amplicon_1_button.config(command=lambda: self.on_button(self.amplicon_1_button,
                                                                     self.amplicon_1_entry, 'amplicon_1_entry'))

        self.amplicon_2_label = ttk.Label(self.master, text='Amplicon2 range: \nformat: XXX-XXX', width=15)
        self.amplicon_2_label.grid(row=12, column=0, rowspan=2)
        self.amplicon_2_entry = ttk.Entry(self.master, width=12)
        self.amplicon_2_entry.grid(row=12, column=1, rowspan=2)
        self.amplicon_2_button = ttk.Button(self.master, text = 'confirm')
        self.amplicon_2_button.grid(row=12, column=2, rowspan=2)
        self.amplicon_2_button.config(command=lambda: self.on_button(self.amplicon_2_button,
                                                                    self.amplicon_2_entry, 'amplicon_2_entry'))

        self.amplicon_3_label = ttk.Label(self.master, text='Amplicon3 range:\nformat: XXX-XXX', width=15)
        self.amplicon_3_label.grid(row=14, column=0, rowspan=2)
        self.amplicon_3_entry = ttk.Entry(self.master, width=12)
        self.amplicon_3_entry.grid(row=14, column=1, rowspan=2)
        self.amplicon_3_button = ttk.Button(self.master, text = 'confirm')
        self.amplicon_3_button.grid(row=14, column=2, rowspan=2)
        self.amplicon_3_button.config(command=lambda: self.on_button(self.amplicon_3_button,
                                                                     self.amplicon_3_entry, 'amplicon_3_entry'))

        self.amplicon_4_label = ttk.Label(self.master, text='Amplicon4 range:\nformat: XXX-XXX', width=15)
        self.amplicon_4_label.grid(row=16, column=0, rowspan=2)
        self.amplicon_4_entry = ttk.Entry(self.master, width=12)
        self.amplicon_4_entry.grid(row=16, column=1, rowspan=2)
        self.amplicon_4_button = ttk.Button(self.master, text = 'confirm')
        self.amplicon_4_button.grid(row=16, column=2, rowspan=2)
        self.amplicon_4_button.config(command=lambda: self.on_button(self.amplicon_4_button,
                                                                     self.amplicon_4_entry, 'amplicon_4_entry'))

        self.amplicon_5_label = ttk.Label(self.master, text='Amplicon5 range:\nformat: XXX-XXX', width=15)
        self.amplicon_5_label.grid(row=18, column=0, rowspan=2)
        self.amplicon_5_entry = ttk.Entry(self.master, width=12)
        self.amplicon_5_entry.grid(row=18, column=1, rowspan=2)
        self.amplicon_5_button = ttk.Button(self.master, text = 'confirm')
        self.amplicon_5_button.grid(row=18, column=2, rowspan=2)
        self.amplicon_5_button.config(command=lambda: self.on_button(self.amplicon_5_button,
                                                                     self.amplicon_5_entry, 'amplicon_5_entry'))

        self.mutation_label = ttk.Label(self.master, text='Target sites in WT \nfrom each amplicon:\nformat: (1) (61 XXX)', width=15)
        self.mutation_label.grid(row=20, column=0, rowspan=3)
        self.mutation_entry = ttk.Entry(self.master, width=29)
        self.mutation_entry.grid(row=20, column=1, rowspan=3, columnspan=2)
        self.mutation_button = ttk.Button(self.master, text = 'confirm')
        self.mutation_button.grid(row=20, column=3, rowspan=3)
        self.mutation_button.config(command=lambda: self.on_button(self.mutation_button, self.mutation_entry, 'mutation_entry'))

        self.submit_button = ttk.Button(self.master, text='      Submit', width=10)
        self.submit_button.grid(row=23, column=1, rowspan=1, columnspan=2)
        self.submit_button.config(command=self.on_submit)

    def openfile(self, button, buttonname):
        filename = askopenfilename(title="Choose file:")
        button.config(text="Checked!")
        button.config(state='disabled')
        self.Dict[buttonname] = filename

    def on_button(self, button, entry, entryname):
        self.Dict[entryname] = entry.get()
        button.config(text='confirmed!', state='disabled')

    def on_submit(self):
        #self.submit_button.config(text='Submitted! Windows will close when finished ...', width=35)
        self.submit_button.config(state='disabled')
        self.master.destroy()

class GuiParam():
    def __init__(self, guiDict):
        self.parseDict(guiDict)
    def parseDict(self, guiDict):
        self.wtfile = None
        self.ngs_data_local = []
        self.condition = []
        self.amplicon  = []
        self.mut_list  = []
        self.mut_pos   = []

def main():
    root = tk.Tk()
    app  = App(root)
    root.mainloop()

    try:
        config = GuiParam(app.Dict)
        config.wtfile = Path(app.Dict['wt_file_button'])

        for item in ['exp_1', 'exp_2', 'exp_0']:
            for key in app.Dict:
                if key.startswith(item):
                    if 'entry' in key:
                        config.condition.append(app.Dict[key])
                    elif 'button' in key:
                        config.ngs_data_local.append(app.Dict[key])
        for item in ['amplicon_1_entry', 'amplicon_2_entry', 'amplicon_3_entry', 'amplicon_4_entry', 'amplicon_5_entry']:
            for key in app.Dict:
                if key == item:
                    config.amplicon.append(tuple(map(int,app.Dict[key].split('-'))))
        if not 'mutation_entry' in app.Dict:
            print("Not enough info provided to GUI, re-try please.")
        else:
            try:
                start = app.Dict['mutation_entry']
                _mut_list = [[y for y in x.split("(") if not y in ['', ' ']] for x in start.split(")") if not x in ['', ' ']]
                _mut_list = [x[0] for x in _mut_list]
                _mut_tuple = [(x+',').replace(' ', ',') for x in _mut_list]
                _mut_tuple = list(map(make_tuple, _mut_tuple))
                mut_tuple = sorted(_mut_tuple, key=lambda item: item[0])
                config.mut_pos = mut_tuple
                config.mut_list = [item for subtuple in mut_tuple for item in subtuple]

                workdir    = Path(Path.cwd()).parents[0]
                param = config
                Path(workdir.joinpath('TemFolder')).mkdir(parents=True, exist_ok=True) # create a temperate folder to contain tem files.

                wrapper.func_wrapper(param, workdir)
            except:
                print("Incorrect information provided to GUI, re-try please.")
    except:
        print("Input format(s) incorrect! Please follow the input format examples and re-try.")


if __name__ == "__main__":
    main()
