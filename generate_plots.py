from bokeh.plotting import figure
from bokeh.charts import Bar
from bokeh.embed import components
import json

doc_head = '<!DOCTYPE html>\n<html lang="en">\n    <head>\n        <meta charset="utf-8">\n        <title>Bokeh Plot</title>\n        \n<link rel="stylesheet" href="https://cdn.pydata.org/bokeh/release/bokeh-0.12.2.min.css" type="text/css" />\n        \n<script type="text/javascript" src="https://cdn.pydata.org/bokeh/release/bokeh-0.12.2.min.js"></script>\n<script type="text/javascript">\n    Bokeh.set_log_level("info");\n</script>\n        <style>\n          html {\n            width: 100%;\n            height: 100%;\n          }\n          body {\n            width: 90%;\n            height: 100%;\n            margin: auto;\n          }\n        </style>\n    </head>\n    <body>\n'
doc_tail = '</body>\n</html>'


bokeh_html = ""

dat=json.loads(open("plot_data.json").readline())


### Sample Missing
fig = figure(title="Sample Missingness", x_axis_label='sorted samples index', y_axis_label='% SNP positions missing')
y = sorted([ x["val"] for x in dat["sampleStats"]["na"]])
x = range(len(y))
fig.line(x,y)   
script, div = components(fig)
bokeh_html += div+script

### Sample Missing
fig = figure(title="Sample Mixedness", x_axis_label='sorted samples index', y_axis_label='% SNP positions mixed')
y = sorted([ x["val"] for x in dat["sampleStats"]["mx"]])
x = range(len(y))
fig.line(x,y)   
script, div = components(fig)
bokeh_html += div+script

### variant Missingness
fig = figure(title="Variant Missingness", x_axis_label='sorted samples index', y_axis_label='% Samples missing SNP data')
y = sorted(dat["varStats"]["na"])
x = range(len(y))
fig.line(x,y)   
script, div = components(fig)
bokeh_html += div+script

### variant Missingness
fig = figure(title="Variant Mixedness", x_axis_label='sorted samples index', y_axis_label='% Samples mixed SNP data')
y = sorted(dat["varStats"]["mx"])
x = range(len(y))
fig.line(x,y)   
script, div = components(fig)
bokeh_html += div+script

### Sample Filtering
bar = Bar(dat["sampleFilter"],values="values",label="labels")
script, div = components(bar)
bokeh_html += div+script

### Variant Filtering
bar = Bar(dat["varFilter"],values="values",label="labels")
script, div = components(bar)
bokeh_html += div+script



open("index.html","w").write(doc_head+"\n"+bokeh_html+"\n"+doc_tail)
