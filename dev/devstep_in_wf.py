"""Dev step downstream of a DataStep in an existing workflow: relation = that step's computed relation.
   STUDIO_TOK=<tokenfile> python dev/devstep_in_wf.py <wf> <upstreamStep> <name> <rows> <cols> <y> <ns> [--filter f=v ...] [prop=value ...]
   rows/cols: comma lists of name[:type]. Prints STEP=<id> (the workflow id is unchanged)."""
import os, sys, uuid, json
import tercen.model.impl as m
from tercen.client.factory import TercenClient
tok = open(os.environ["STUDIO_TOK"]).read().strip()
c = TercenClient("http://127.0.0.1:5402"); c.userService.tercenClient.token = tok; c.httpClient.authorization = tok
wf_id, up_id, name, rows, cols, y, ns = sys.argv[1:8]; rest = sys.argv[8:]
filters = []; props = {}
i = 0
while i < len(rest):
    if rest[i] == "--filter": filters.append(rest[i + 1].split("=", 1)); i += 2
    else: k, v = rest[i].split("=", 1); props[k] = v; i += 1
wf = c.workflowService.get(wf_id); up = next(s for s in wf.steps if s.id == up_id)
task = c.taskService.get(up.state.taskId); rel = task.computedRelation
def rect(x, y, w=200.0, h=55.0):
    r = m.Rectangle(); r.topLeft = m.Point(); r.topLeft.x = x; r.topLeft.y = y; r.extent = m.Point(); r.extent.x = w; r.extent.y = h; return r
def typed(spec):
    out = []
    for x in spec.split(","):
        if not x: continue
        n, _, t = x.partition(":"); out.append((n, t or "string"))
    return out
def factor(n, t):
    f = m.Factor(); f.name = n; f.type = t; return f
def gf(n, t):
    g = m.GraphicalFactor(); g.factor = factor(n, t); g.rectangle = rect(0.0, 0.0, 60.0, 30.0); return g
def ctable(fs):
    t = m.CrosstabTable(); t.cellSize = 250.0; t.offset = 0; t.nRows = 0; t.graphicalFactors = fs; t.rectangleSelections = []; return t
def port(cls, nm):
    p = cls(); p.id = str(uuid.uuid4()); p.name = nm; p.linkType = "relation"; return p
ds = m.DataStep(); ds.id = str(uuid.uuid4()); ds.name = name; ds.groupId = ""; ds.description = ""; ds.parentDataStepId = ""
ip = port(m.InputPort, "data"); ds.inputs = [ip]; ds.outputs = [port(m.OutputPort, "data")]; ds.rectangle = rect(up.rectangle.topLeft.x + 300.0, up.rectangle.topLeft.y + 150.0)
ds.state = m.StepState(); ds.state.taskId = ""; ds.state.taskState = m.InitState()
ct = m.Crosstab(); ct.taskId = ""; ct.axis = m.XYAxisList(); ct.axis.rectangleSelections = []; ct.axis.xyAxis = []
ct.columnTable = ctable([gf(n, t) for n, t in typed(cols)] or [gf("", "string")]); ct.rowTable = ctable([gf(n, t) for n, t in typed(rows)] or [gf("", "string")])
ct.filters = m.Filters(); ct.filters.removeNaN = False; ct.filters.namedFilters = []
if filters:
    nf = m.NamedFilter(); nf.name = "dev"; nf.logical = "and"; nf.meta = []; nf.isNot = False; nf.filterExprs = []
    for fn, fv in filters:
        fe = m.FilterExpr(); fe.factor = factor(fn, "string"); fe.filterOp = "equals"; fe.stringValue = fv; fe.preProcessors = []; nf.filterExprs.append(fe)
    ct.filters.namedFilters = [nf]
st = m.OperatorSettings(); st.namespace = ns; st.environment = []
ref = m.OperatorRef(); ref.name = "dev"; ref.version = "dev"; ref.operatorId = ""; ref.operatorKind = ""; ref.url = m.Url(); ref.url.uri = ""; ref.propertyValues = []
for k, v in props.items():
    pv = m.PropertyValue(); pv.name = k; pv.value = v; ref.propertyValues.append(pv)
ref.operatorSpec = m.OperatorSpec(); ref.operatorSpec.inputSpecs = []; ref.operatorSpec.outputSpecs = []
st.operatorRef = ref; ct.operatorSettings = st; ds.model = ct
link = m.Link(); link.id = str(uuid.uuid4()); link.inputId = ip.id; link.outputId = up.outputs[0].id
wf.steps = list(wf.steps) + [ds]; wf.links = list(wf.links) + [link]; c.workflowService.update(wf)
q = m.CubeQuery(); q.relation = rel; q.colColumns = [factor(n, t) for n, t in typed(cols)]; q.rowColumns = [factor(n, t) for n, t in typed(rows)]
aq = m.CubeAxisQuery(); aq.chartType = "point"; aq.pointSize = 4
aq.xAxis = factor("", "string"); aq.yAxis = factor(y, "double"); aq.colors = []; aq.errors = []; aq.labels = []; aq.preprocessors = []
aq.xAxisSettings = m.AxisSettings(); aq.xAxisSettings.meta = []; aq.yAxisSettings = m.AxisSettings(); aq.yAxisSettings.meta = []
q.axisQueries = [aq]; q.filters = ct.filters; q.operatorSettings = st
t2 = m.CubeQueryTask(); t2.state = m.InitState(); t2.owner = wf.acl.owner; t2.projectId = wf.projectId; t2.query = q
t2 = c.taskService.create(t2); c.taskService.runTask(t2.id); t2 = c.taskService.waitDone(t2.id)
print("CUBE", type(t2.state).__name__, getattr(t2.state, "reason", ""))
wf = c.workflowService.get(wf_id); ds = next(s for s in wf.steps if s.id == ds.id); ds.model.taskId = t2.id
axis = json.loads('{"kind":"XYAxisList","xyAxis":[{"kind":"XYAxis","chart":{"kind":"ChartPoint","name":"","pointSize":4,"properties":{"kind":"Properties","properties":[],"propertyValues":[]}},"xAxis":{"kind":"Axis","axisExtent":{"x":80.0,"y":30.0,"kind":"Point"},"axisSettings":{"kind":"AxisSettings","meta":[]},"graphicalFactor":{"kind":"GraphicalFactor","factor":{"kind":"Factor","name":"","type":"string"},"rectangle":{"kind":"Rectangle","extent":{"x":0.0,"y":0.0,"kind":"Point"},"topLeft":{"x":0.0,"y":0.0,"kind":"Point"}}}},"yAxis":{"kind":"Axis","axisExtent":{"x":80.0,"y":30.0,"kind":"Point"},"axisSettings":{"kind":"AxisSettings","meta":[]},"graphicalFactor":{"kind":"GraphicalFactor","factor":{"kind":"Factor","name":"value","type":"double"},"rectangle":{"kind":"Rectangle","extent":{"x":0.0,"y":0.0,"kind":"Point"},"topLeft":{"x":0.0,"y":0.0,"kind":"Point"}}}},"colors":{"kind":"Colors","factors":[],"palette":{"kind":"CategoryPalette","backcolor":0,"colorList":{"kind":"ColorList","name":""},"properties":[],"stringColorElements":[]}},"errors":{"kind":"Errors","factors":[]},"labels":{"kind":"Labels","factors":[]},"taskId":"","preprocessors":[]}],"rectangleSelections":[]}')
axis["xyAxis"][0]["yAxis"]["graphicalFactor"]["factor"]["name"] = y; axis["xyAxis"][0]["taskId"] = t2.id
ds.model.axis = m.XYAxisList(axis); c.workflowService.update(wf)
print(f"STEP={ds.id}")
