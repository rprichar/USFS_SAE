from arcgis.features import FeatureLayer
from arcgis.geometry import Geometry, MultiPoint
from arcgis.gis import GIS
from arcgis.raster import ImageryLayer
import arcpy
import numpy as np
import os
import pandas as pd
import random
import time
arcpy.env.overwriteOutput = True

# https://pro.arcgis.com/en/pro-app/latest/arcpy/geoprocessing_and_python/a-template-for-python-toolboxes.htm

class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Toolbox"
        self.alias = ""

    # List of tool classes associated with this toolbox
        self.tools = [GetSamples]


class GetSamples(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "GetSamples"
        self.description = "GetSamples"
        self.canRunInBackground = True

    def getParameterInfo(self):
        param_input_aoi = arcpy.Parameter(
            displayName="Input Area of Interest (polygon)",
            name="input_aoi",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Input"
        )

        param_field = arcpy.Parameter(
            displayName="Area of Interest Field",
            name="field",
            datatype="String",
            parameterType="Required",
            direction="Input"
        )

        param_image_service = arcpy.Parameter(
            displayName="Image Service URL",
            name="image_service",
            datatype="DEImageServer",
            parameterType="Required",
            direction="Input"
        )

        param_k_no = arcpy.Parameter(
            displayName="Maximum Number of Neighbors (int)",
            name="k_no",
            datatype="GPLong",
            parameterType="Required",
            direction="Input"
        )

        param_sampling_factor = arcpy.Parameter(
            displayName="Proportion of Sample Pixels (0-1)",
            name="sampling_factor",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input"
        )

        param_output_csv_dir = arcpy.Parameter(
            displayName="Output directory for the .csv file",
            name="output_csv_dir",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input"
        )

        param_output_csv = arcpy.Parameter(
            displayName="Full output path and name of the .csv file",
            name="output_csv",
            datatype="String",
            parameterType="Derived",
            direction="Output"
        )

        params = [param_input_aoi,
                param_field,
                param_image_service,
                param_k_no,
                param_sampling_factor,
                param_output_csv_dir,
                param_output_csv
                ]
        return params

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
            parameters[1].filter.list = [f.name for f in arcpy.ListFields(parameters[0].value)]
            return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        # Read in the parameters
        param_input_aoi = parameters[0].valueAsText
        param_field = parameters[1].valueAsText
        param_image_service = parameters[2].valueAsText
        param_k_no = parameters[3].valueAsText
        param_sampling_factor = parameters[4].valueAsText
        param_output_csv_dir = parameters[5].valueAsText

        gis = GIS("pro")
        img_lyr = ImageryLayer(param_image_service,  gis=gis)

        class LocalTimer(object):
            def __init__(self):
                # Report start time
                self.starttime = time.time()
                print("Start Time: {0}".format(self.formattime(self.starttime)))

            def formattime(self,ltime):
                return time.asctime(time.localtime(ltime))

            def reportelapsedtime(self,etime):
                stime=int(etime%60)
                if (stime == 1):
                    sstring = " 1 second"
                else:
                    sstring = " {0} seconds".format(stime)

                if (etime < 60):
                    mstring = ""
                    hstring = ""
                else:
                    mtime = int((etime%3600)/60)
                    if (mtime == 1):
                        mstring = " 1 minute"
                    else:
                        mstring = " {0} minutes".format(mtime)

                    if (etime <3600):
                        hstring = ""
                    else:
                        htime = int(etime/3600)
                        if (htime == 1):
                            hstring = " 1 hour"
                        else:
                            hstring = " {0} hours".format(htime)

                return "{0}{1}{2}".format(hstring, mstring, sstring)

            def delete(self):
                # Report time taken
                if self.starttime > 0:
                    endtime = time.time()
                    elapsedtime = endtime - self.starttime
                    print("End Time: {0} (Elapsed Time:{1})".format(self.formattime(endtime), self.reportelapsedtime(elapsedtime)))
                    self.starttime = 0

            def __enter__(self):
                return self

            def __exit__(self,type,value,traceback):
                self.delete()

            def __del__(self):
                self.delete()

        class arcpy_gp():
            def __init__(self, aoi, sampling_factor, group_by_field):
                self.aoi = aoi
                self.sampling_factor = sampling_factor
                self.group_by_field = group_by_field
                
            def add_fields(self, in_feature, field_desc):
                arcpy.management.AddFields(in_feature, field_desc)

            def apply_sampling_factor(self, row_count, sampling_factor):
                samples = int(row_count*sampling_factor)
                return(samples)
                
            def calculate_field(self, in_feature, field, expression):
                arcpy.management.CalculateField(in_feature, field, expression)
                
            def clip(self, in_feature, clip_feature, out_feature):
                arcpy.analysis.Clip(in_feature, clip_feature, out_feature)

            def create_fishnet_30mgrid(self, in_feature, origin_coordinate, y_axis_coordinate, opposite_coordinate, full_extent):
                arcpy.management.CreateFishnet(in_feature, 
                                    origin_coordinate, 
                                    y_axis_coordinate, 
                                    30, 30, None, None, 
                                    opposite_coordinate, 
                                    "LABELS", 
                                    full_extent, 
                                    "POLYGON")
                
            def create_sel_statement(self, row_count, samples):
                randomlist = [str(i) for i in random.sample(range(0, row_count), samples)]
                randomlist = ['OBJECTID = ' + i for i in randomlist]
                sel_statement = " OR ".join(randomlist)
                return(sel_statement)
            
            def delete(self, in_feature):
                arcpy.management.Delete(in_feature)
                
            def define_projection_WGS84(self, in_feature):
                arcpy.management.DefineProjection(in_feature, 
                                            'PROJCS["WGS_1984_Web_Mercator_Auxiliary_Sphere",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Mercator_Auxiliary_Sphere"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",0.0],PARAMETER["Standard_Parallel_1",0.0],PARAMETER["Auxiliary_Sphere_Type",0.0],UNIT["Meter",1.0]]')

            def get_row_count(self, in_feature):
                row_count = arcpy.management.GetCount(in_feature)
                row_count = int(row_count[0])
                return(row_count)
            
            def list_fields(self,):
                field_names = [f.name for f in arcpy.ListFields(self.aoi)]
                return(field_names)    
            
            def select_and_copy(self, in_feature, workspace, out_feature_name, sel_statement):
                sel = arcpy.management.SelectLayerByAttribute(in_feature, "NEW_SELECTION", sel_statement, None)
                arcpy.conversion.FeatureClassToFeatureClass(sel, workspace, out_feature_name) 
                
            def sort_peano_curve(self, in_feature, out_feature):
                arcpy.Sort_management(in_feature, out_feature, [["Shape", "ASCENDING"]], "PEANO")
                

            def sae_execute_preprocessing(self,):
                for row in arcpy.da.SearchCursor(self.aoi, ["SHAPE@",self.group_by_field]):
                    polygon_id = str(row[1])
                    polygon_id = ''.join(i for i in polygon_id if i.isalnum())
                    fishnet = self.aoi + '_' + polygon_id + '_fishnet'
                    fishnet_label = self.aoi + '_' + polygon_id + '_fishnet_label' 
                    fishnet_label_aoi = self.aoi + '_' + polygon_id + '_fishnet_label_aoi'
                    fishnet_selection = self.aoi + '_' + polygon_id + '_fishnet_selection'
                    sampling_peano = self.aoi + '_' + polygon_id + '_sampling_peano'

                    extent = row[0].extent
                    origin_coordinate =  '{} {}'.format(extent.XMin, extent.YMin)
                    y_axis_coordinate =  '{} {}'.format(extent.XMin, extent.YMin +10)
                    opposite_coordinate =  '{} {}'.format(extent.XMax, extent.YMax)
                    full_extent = '{} {} {} {}'.format(extent.XMin, extent.YMin, extent.XMax, extent.YMax)

                    self.create_fishnet_30mgrid(fishnet, origin_coordinate, y_axis_coordinate, opposite_coordinate, full_extent)
                    self.define_projection_WGS84(fishnet_label)
                    self.clip(fishnet_label, self.aoi, fishnet_label_aoi)
                    row_count = self.get_row_count(fishnet_label_aoi)
                    samples = self.apply_sampling_factor(row_count, self.sampling_factor)
                    sel_statement = self.create_sel_statement(row_count, samples)
                    #arcpy.AddMessage("sel_statement:{0}".format(sel_statement))
                    self.select_and_copy(fishnet_label_aoi, os.path.dirname(fishnet_selection), os.path.basename(fishnet_selection), sel_statement)
                    self.sort_peano_curve(fishnet_selection, sampling_peano)
                    self.calculate_field(sampling_peano, "polygon_id", '"' + polygon_id + '"')
                    self.delete(fishnet)
#                    self.delete(fishnet_label)
#                    self.delete(fishnet_label_aoi)
#                    self.delete(fishnet_selection)
                    
            def walk_workspace_return_feature_class_with_selection(self):
                feature_classes = []
                walk = arcpy.da.Walk(os.path.dirname(self.aoi), datatype="FeatureClass", type="Point")
                for dirpath, dirnames, filenames in walk:
                    for filename in filenames:
                        if '_sampling_peano' in filename:
                            feature_classes.append(os.path.join(dirpath, filename))
                return(feature_classes)

        class df():
            def __init__(self, input_list, polygon_id, k_no, output_csv):
                self.input_list = input_list
                self.polygon_id = polygon_id
                self.k_no = k_no
                self.output_csv = output_csv

            def create_from_list(self, **kwargs):
                df = pd.DataFrame(data=self.input_list, **kwargs)
                return(df)

            def calculate_field(self, df, field_name, field_value):
                df[field_name] = field_value
                return(df)
                
            def drop_columns_on_k(self, df, k_no):
                keep_list = ['x', 'y']
                for i in range(2, k_no+2):
                    keep_list.append('k_'+ str(i))
                df = df[df.columns.intersection(keep_list)]
                return(df)
            
            def drop_NA(self, df):
                df = df.dropna()
                return(df)

            def label_xyk(self, row_count):
                label_list = ['x', 'y']
                for i in range(1, row_count-1):
                    label_list.append('k_'+ str(i))
                return(label_list)    

            def return_row_count(self, df):
                row_count = df.shape[1]
                return(row_count)
            
            def return_col_count(self, df):
                col_count = df.shape[0]
                return(col_count)
            
            def save_to_csv(self,df):
                df.to_csv(self.output_csv, index=False)
            
            def execute_get_samples_to_csv(self):
                df = self.create_from_list()
                row_count = self.return_row_count(df)
                label_list = self.label_xyk(row_count)
                df=self.create_from_list(columns=label_list)
                df=self.drop_columns_on_k(df, self.k_no)
                col_count_full = self.return_col_count(df)
                df=self.calculate_field(df,'area_full_m2', col_count_full*900)
                df=self.drop_NA(df)
                col_count_valid = self.return_col_count(df)
                df=self.calculate_field(df,'area_valid_m2', col_count_valid*900)
                df=self.calculate_field(df,'area_null_m2', (col_count_full-col_count_valid)*900)
                df=self.calculate_field(df,'polygon_id',self.polygon_id)
                self.save_to_csv(df)
                return(df)        

        def call_url(img_lyr, mp, attempts=5, chunk=None):    
            ret_list = []
            pixels   = []
            geom_len = len(mp['points'])
            for attempt in range(1, attempts+1):
                new_ret = []
                try:
                    new_ret = img_lyr.get_samples(mp)
                    if len(new_ret) < geom_len:
                        if len(new_ret) > len(pixels):
                            pixels = new_ret

                        if attempt == attempts:
                            print ('Failed to get all values after {} attempts'.format(attempts))
                        else:
                            print ('Not enough values found {} out of {} after {} attempts, retrying...'.format(len(new_ret), geom_len, attempt))
                    else:
                        pixels = new_ret

                        if attempt > 1:
                            print ('All values found after {} attempts'.format(attempt))
                        break
                except:
                    if attempt == attempts:
                        if len(pixels) == 0:
                            print ('Error: failed to get any values after {} attempts'.format(attempts))
                            return None
                    else:
                        print ('Error: failed to get values on attempt {}, retrying...'.format(attempt))
                    pass

            try:
                for pixel in pixels:
                    px      = pixel['location']['x']
                    py      = pixel['location']['y']
                    pvalstr = pixel['value']
                    pvals   = [float(s) for s in pvalstr.split()]
                    pline = [px, py]
                    pline.extend(pvals)
                    ret_list.append(pline)
            
                plen = len(pixels)
                pmessage = 'Found {} values for {} points'.format(plen,geom_len)
                if chunk:
                    pmessage = pmessage + ' on batch {}'.format(chunk)
                if plen < geom_len:
                    pmessage = pmessage + ' - ERROR: {} missing after {} attempts'.format(geom_len-plen,attempt)
                print(pmessage)

                return ret_list
            
            except:
                print('Error: cannot compile results')
                return None

        param_k_no = int(param_k_no)
        #arcpy.AddMessage("param_sampling_factor: {0}".format(param_sampling_factor))
        sampling_factor = float(param_sampling_factor)
        #arcpy.AddMessage("sampling_factor: {0}".format(sampling_factor))

        polygon_id = param_field
        #arcpy.AddMessage("polygon_id: {0}".format(polygon_id))
        arcpy_gp = arcpy_gp(param_input_aoi, sampling_factor, polygon_id)
        arcpy_gp.sae_execute_preprocessing()
        arcpy.AddMessage("Preprocessing of input data complete")
        sampling_list = arcpy_gp.walk_workspace_return_feature_class_with_selection()
        arcpy.AddMessage("Submitting the following sampling files:{0}".format(sampling_list))

        for sample in sampling_list:
            with arcpy.da.SearchCursor(sample, ["SHAPE@XY",'polygon_id']) as cur:
                points = []
                for row in cur:
                    x, y = row[0]
                    points.append([x, y])
                    polygon_id = str(row[1])

            # Reproject to service projection system
            in_geom = Geometry({'points':points, 'spatialReference':{'wkid':3857}})
            #print('Number of records = {}'.format(len(in_geom['points'])))
            arcpy.AddMessage('Number of records = {}'.format(len(in_geom['points'])))

            ##################################################################################

            # Set batch size
            bsize   = 200

            mp_len  = len(in_geom['points'])
            mp_list = []

            for olen in range(0, mp_len, bsize):
                nlen = olen + bsize
                if nlen > mp_len:
                    nlen = mp_len

                geom = {'points' :in_geom['points'][olen:nlen], "spatialReference": {
                    "wkid": 102100,
                    "latestWkid": 3857
                }}
                mp_list.append(MultiPoint(geom))

            # How many chunks do we now have
            #print('Number of calls: {}'.format(len(mp_list)))
            arcpy.AddMessage('Number of calls: {}'.format(len(mp_list)))

            # Variables to control execution
            #img_lyr  = img_lyr
            attempts = 1
            nthreads = 4

            # Output list for writing to CSV
            plist    = []

            # Execute calls
            with LocalTimer() as localtimer:
                olists   = []

                if nthreads > 1:
                    from concurrent import futures

                    #print ('Executing calls in {} parallel threads'.format(nthreads))
                    arcpy.AddMessage('Executing calls in {} parallel threads'.format(nthreads))
                    with futures.ThreadPoolExecutor(max_workers=nthreads) as executor:
                        calls = {executor.submit(call_url, img_lyr, mp, attempts): mp for mp in mp_list}
                        for call in futures.as_completed(calls):
                            try:
                                olist = call.result()
                            except Exception as exc:
                                #print('Exception generated: %s' % (exc))
                                arcpy.AddMessage('Exception generated: %s' % (exc))
                            else:
                                olists.append(olist)
                else:
                    #print ('Executing calls in a single thread')
                    arcpy.AddMessage('Executing calls in a single thread')
                    olists = [call_url(param_image_service, mp_list[mp] , attempts, mp+1) for mp in range(len(mp_list))]

                #print('{} values returned'.format(len(olists)))
                arcpy.AddMessage('{} values returned'.format(len(olists)))

                # Combine lists
                for olist in olists:
                    if olist is not None:
                        plist.extend(olist)

                #print('{} data values found - {} missing'.format(len(plist),len(in_geom['points'])-len(plist)))
                arcpy.AddMessage('{} data values found - {} missing'.format(len(plist),len(in_geom['points'])-len(plist)))

            ##################################################################################
                param_output_csv = os.path.join(param_output_csv_dir, os.path.basename(sample)) + '.csv'
                df_inst = df(plist, polygon_id, param_k_no, param_output_csv)
                df_inst.execute_get_samples_to_csv()
                arcpy.AddMessage('{} file written'.format(param_output_csv))
#            arcpy_gp.delete(sample)         