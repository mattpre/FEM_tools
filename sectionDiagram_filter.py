from vtkmodules.vtkCommonDataModel import vtkDataSet, vtkPolyData, vtkPlane, vtkVector3d, vtkCellArray, vtkQuad, vtkHexahedron
from vtkmodules.vtkCommonCore import vtkFloatArray, vtkDataArraySelection
from vtkmodules.vtkFiltersCore import vtkCutter
from vtkmodules.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtkmodules.numpy_interface import dataset_adapter as dsa

from paraview.util.vtkAlgorithm import smproxy, smproperty, smdomain

def GetQuad(p0,p1,ids,n,t,val,points):
    p2 = [p1[k]+n[k]*val for k in range(3)]
    p3 = [p0[k]+n[k]*val for k in range(3)]
    id2 = points.InsertNextPoint(p2)
    id3 = points.InsertNextPoint(p3)
    quad = vtkQuad()
    quad.GetPointIds().SetId(0,ids.GetId(0))
    quad.GetPointIds().SetId(1,ids.GetId(1))
    quad.GetPointIds().SetId(2,id2)
    quad.GetPointIds().SetId(3,id3)

    return quad

def GetQuads(p0,p1,n,n2,t,val,points):
    p00 = [p0[k]-n2[k] for k in range(3)]
    p01 = [p0[k]+n2[k] for k in range(3)]
    p10 = [p1[k]-n2[k] for k in range(3)]
    p11 = [p1[k]+n2[k] for k in range(3)]
    p20 = [p0[k]+n[k]*val-n2[k] for k in range(3)]
    p21 = [p0[k]+n[k]*val+n2[k] for k in range(3)]
    p30 = [p1[k]+n[k]*val-n2[k] for k in range(3)]
    p31 = [p1[k]+n[k]*val+n2[k] for k in range(3)]
    id0 = points.InsertNextPoint(p00)
    id1 = points.InsertNextPoint(p01)
    id2 = points.InsertNextPoint(p10)
    id3 = points.InsertNextPoint(p11)
    id4 = points.InsertNextPoint(p20)
    id5 = points.InsertNextPoint(p21)
    id6 = points.InsertNextPoint(p30)
    id7 = points.InsertNextPoint(p31)
    q0 = vtkQuad()
    q0.GetPointIds().SetId(0,id0)
    q0.GetPointIds().SetId(1,id1)
    q0.GetPointIds().SetId(2,id3)
    q0.GetPointIds().SetId(3,id2)
    q1 = vtkQuad()
    q1.GetPointIds().SetId(0,id4)
    q1.GetPointIds().SetId(1,id5)
    q1.GetPointIds().SetId(2,id7)
    q1.GetPointIds().SetId(3,id6)
    q2 = vtkQuad()
    q2.GetPointIds().SetId(0,id2)
    q2.GetPointIds().SetId(1,id3)
    q2.GetPointIds().SetId(2,id7)
    q2.GetPointIds().SetId(3,id6)
    q3 = vtkQuad()
    q3.GetPointIds().SetId(0,id0)
    q3.GetPointIds().SetId(1,id1)
    q3.GetPointIds().SetId(2,id5)
    q3.GetPointIds().SetId(3,id4)
    q4 = vtkQuad()
    q4.GetPointIds().SetId(0,id0)
    q4.GetPointIds().SetId(1,id4)
    q4.GetPointIds().SetId(2,id6)
    q4.GetPointIds().SetId(3,id2)
    q5 = vtkQuad()
    q5.GetPointIds().SetId(0,id1)
    q5.GetPointIds().SetId(1,id5)
    q5.GetPointIds().SetId(2,id7)
    q5.GetPointIds().SetId(3,id3)
    

    return [q0,q1,q2,q3,q4,q5]

@smproxy.filter(label="Section Diagram Filter")
##@smproperty.input(name="Input")
@smproperty.xml("""
        <InputProperty command="SetInputConnection"
                      name="Input">
         <ProxyGroupDomain name="groups">
           <Group name="sources" />
           <Group name="filters" />
         </ProxyGroupDomain>
         <DataTypeDomain name="input_type">
           <DataType value="vtkDataSet" />
         </DataTypeDomain>
         <InputArrayDomain attribute_type="cell"
                           name="input_array1"
                           number_of_components="1"
                           optional="1" />
         <Documentation>This property specifies the input to the Glyph
 filter.
         This is the dataset from which the points are selecetd to be
 glyphed.
         </Documentation>
       </InputProperty>
""")
class SectionDiagramFilter(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self)

        self.normal = vtkVector3d(1,0,0)
        self.origin = [0,0,0]
        self.scale = 1
        self.component = 0
        self.selected_array = None
        self.thickness = 0

        self.activearray = 0

        self._arrayselection = vtkDataArraySelection()

    def RequestData(self, request, inInfo, outInfo):
        input0 = dsa.WrapDataObject(vtkDataSet.GetData(inInfo[0]))

        cd = input0.VTKObject.GetCellData()
        for kcd in range(cd.GetNumberOfArrays()):
            self._arrayselection.AddArray(cd.GetArray(kcd).GetName())

        output = dsa.WrapDataObject(vtkDataSet.GetData(outInfo))

        plane = vtkPlane()
        plane.SetOrigin(self.origin)
        plane.SetNormal(self.normal)

        cut = vtkCutter()
        cut.SetInputData(input0.VTKObject)
        cut.SetCutFunction(plane)
        cut.Update()

        cut_out = cut.GetOutput()
        if (cut_out.GetNumberOfCells()):
            points = cut_out.GetPoints()
            cells = vtkCellArray()
            array = vtkFloatArray()
            array.SetNumberOfComponents(1)

            array.SetName(self.selected_array)        
            arr = cut_out.GetCellData().GetArray(self.selected_array)
            max_val = max([abs(arr.GetTuple(kk)[self.component]) for kk in range(arr.GetNumberOfTuples())])
            sc = 0.02*cut_out.GetLength()/max_val
            for kc in range(cut_out.GetNumberOfCells()):
                cell = cut_out.GetCell(kc)
                ids = cell.GetPointIds()
                p0 = cell.GetPoints().GetPoint(0)
                p1 = cell.GetPoints().GetPoint(1)
                t = vtkVector3d(p1[0]-p0[0],p1[1]-p0[1],p1[2]-p0[2])
                t.Normalize()
                n = t.Cross(self.normal)
                n2 = vtkVector3d(self.normal[0]*self.thickness,
                                 self.normal[1]*self.thickness,
                                 self.normal[2]*self.thickness)
                val = arr.GetTuple(kc)[self.component]
                if self.thickness==0:
                    quad = GetQuad(p0,p1,ids,n,t,val*self.scale*sc,points)
                    cells.InsertNextCell(quad)
                    array.InsertNextTuple1(val)
                else:
                    quads = GetQuads(p0,p1,n,n2,t,val*self.scale*sc,points)
                    for quad in quads:
                        cells.InsertNextCell(quad)
                        array.InsertNextTuple1(val)
            output.SetPolys(cells)
            output.SetPoints(points)
            output.GetCellData().AddArray(array)

        return 1

    @smproperty.doublevector(name="Direction", default_values=[1, 0, 0])
    @smdomain.doublerange()
    def SetDirection(self, x, y, z):
        self.normal = vtkVector3d(x,y,z)
        self.Modified()

    @smproperty.doublevector(name="Origin", default_values=[0, 0, 0])
    @smdomain.doublerange()
    def SetOrigin(self, x, y, z):
        self.origin = [x,y,z]
        self.Modified()

##    @smproperty.doublevector(name="Scale", default_values=1)
##    @smdomain.xml(""" <DoubleRangeDomain min="1" max="100">
##                    </DoubleRangeDomain>""")
##    @smdomain.doublerange(min=0, max=10)
    @smproperty.xml("""
        <DoubleVectorProperty name="Scale"
                              number_of_elements="1"
                              default_values="1"
                              command="SetScale">
            <DoubleRangeDomain min="1" max="100">
            </DoubleRangeDomain>
        </DoubleVectorProperty>
                              """)
    def SetScale(self, x):
        self.scale = x
        self.Modified()

##    @smproperty.stringvector(name="Cell Arrays", number_of_elements="1")
##    @smdomain.xml("""
##        <ArrayListDomain name="array_list"
##                         attribute_type="Scalars"
##                         input_domain_name="inputs_array"
##                         >
##          <RequiredProperties>
##            <Property name="Input"
##                      function="Input" />
##          </RequiredProperties>
##        </ArrayListDomain>
##        <Documentation>
##            Some doc
##        </Documentation>
##        <Hints>
##            <ShowInMenu category="Extensions" />
##        </Hints>
##        """)
    @smproperty.xml("""
        <StringVectorProperty name="Cell Arrays"
                              number_of_elements="1"
                              command="GetDataArraySelection">
            <ArrayListDomain name="array_list"
                         attribute_type="Scalars"
                         input_domain_name="inputs_array"
                         >
              <RequiredProperties>
                  <Property name="Input"
                      function="Input" />
              </RequiredProperties>
            </ArrayListDomain>
            <Documentation>
                Some doc
            </Documentation>
            <Hints>
                <ShowInMenu category="Extensions" />
            </Hints>
        </StringVectorProperty>
        """)
    def GetDataArraySelection(self,name):
        self.selected_array = name
        self.Modified()

    @smproperty.intvector(name="Component", default_values=0)
    @smdomain.intrange(min=0, max=2)
    def SetComponent(self, x):
        self.component = x
        self.Modified()

    @smproperty.doublevector(name="Thickness", default_values=0)
    @smdomain.doublerange(min=0, max=10)
    def SetThickness(self, x):
        self.thickness = x
        self.Modified()
