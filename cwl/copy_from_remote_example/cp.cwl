cwlVersion: v1.0
class: CommandLineTool
baseCommand: cp
inputs:
  src_file:
    type: File
    inputBinding:
      position: 1
  dest_file:
    type: string
    inputBinding:
      position: 2

outputs:
  example_out:
    type: File
    outputBinding:
      glob: $(inputs.dest_file)
