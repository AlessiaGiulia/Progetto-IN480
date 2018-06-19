
using Plots
using Distributions
using Base.Test

#_________________________incremento input con PyCall_____________________
using PyCall
@pyimport larlib as lar
l=[]
for i in range(0,3)
    p=PyObject(lar.larCuboids([10^i,10^i]))
    push!(l,Tuple(append!([map(collect,PyObject(p[1]))],[PyObject(p[2])])))
    p=PyObject(lar.larCuboids([2*10^i,2*10^i]))
    push!(l,Tuple(append!([map(collect,PyObject(p[1]))],[PyObject(p[2])])))
    p=PyObject(lar.larCuboids([5*10^i,5*10^i]))
    push!(l,Tuple(append!([map(collect,PyObject(p[1]))],[PyObject(p[2])])))
end

#____________________________________funzione tempo_______________________________

function Time(f,args)
    @elapsed f(args...)
    t=[]
    for i in range(1,10)
	push!(t,@elapsed f(args...))
    end
    m=mean(t)
    return m
end

#_____________________________incremento input tesla____________________________________

function addn2D(n,model)
  body=[]
  for i in range(1,n)
    el=[]
    matrix=rand(1:3)
    if matrix==1
      x=rand(1:10)/10
      y=rand(1:10)/10
      append!(el,larApply(t(x,y))(model))
      append!(body,[el])
    elseif matrix ==2
      x=rand(1:10)/10
      y=rand(1:10)/10
      append!(el,larApply(s(x,y))(model))
      append!(body,[el])
    elseif matrix ==3
      x=rand(1:10)
      append!(el,larApply(r(pi/x))(model))
      append!(body,[el])
    end
  end
  a=Struct(body)
  return a
end

#_____________________checkStruct_______________________________

function checkStruct(lst)
  obj = lst[1]
  if isa(obj,Matrix)
      dim=size(obj)[1]-1
  elseif(isa(obj,Tuple) || isa(obj,Array))
      dim=length(obj[1][1])
  elseif isa(obj,Struct)
      dim=length(obj.box[1])
  end 
  return dim
end

#____________________pcheckStruct_____________________________

function pcheckStruct(lst)
  obj = lst[1]
  if isa(obj,Matrix)
    dim=size(obj)[1]-1 
  elseif(isa(obj,Tuple) || isa(obj,Array))
    dim=length(obj[1][1])
  elseif isa(obj,pStruct)
    dim=length(obj.box[1])
  end
  return dim
end

#________________"checkStruct tests"______________________________

@testset "checkStruct" begin
  list=([[0.575,-0.175],[0.575,0.175],[0.925,-0.175],[0.925,0.175]],[[0,1,2,3]])
@test checkStruct(list)==length(list[1][1][1])
@test typeof(checkStruct(list))==Int
end

#__________________________larApply______________________________

#rotazione

  function r(args...)
    args = collect(args)
    n = length(args)
    if n == 1 # rotation in 2D
      angle = args[1]; COS = cos(angle); SIN = sin(angle)
      mat = eye(3)
      mat[1,1] = COS;    mat[1,2] = -SIN;
      mat[2,1] = SIN;    mat[2,2] = COS;
    end
    if n == 3 # rotation in 3D
      mat = eye(4)
      angle = norm(args); axis = normalize(args)
      COS = cos(angle); SIN= sin(angle)
      if axis[2]==axis[3]==0.0    # rotation about x
	mat[2,2] = COS;    mat[2,3] = -SIN;
	mat[3,2] = SIN;    mat[3,3] = COS;
      elseif axis[1]==axis[3]==0.0   # rotation about y
	mat[1,1] = COS;    mat[1,3] = SIN;
	mat[3,1] = -SIN;    mat[3,3] = COS;
      elseif axis[1]==axis[2]==0.0    # rotation about z
	mat[1,1] = SIN;    mat[1,2] = -SIN;
	mat[2,1] = COS;    mat[2,2] = COS;
      else
	I=eye(3); u=axis
	Ux=[0 -u[3] u[2] ; u[3] 0 -u[1] ;  -u[2] u[1] 1]
	UU =[u[1]*u[1]    u[1]*u[2]   u[1]*u[3];
	    u[2]*u[1]    u[2]*u[2]   u[2]*u[3];
	    u[3]*u[1]    u[3]*u[2]   u[3]*u[3]]
	mat[1:3,1:3]=COS*I+SIN*Ux+(1.0-COS)*UU
      end
    end
  return mat
  end

#traslazione

 function t(args...)
   d=length(args)
   mat=eye(d+1)
   for k in range(1,d)
     mat[k,d+1]=args[k]
   end
   return mat
 end

#scala

function s(args...)
  d=length(args)
  mat=eye(d+1)
  for k in range(1,d)
    mat[k,k]=args[k]
  end
  return mat
end

#larApply

function larApply(affineMatrix)
  function larApply0(model)
    if length(model)==2
      V,CV=model
    elseif length(model)==3
      V,CV,FV = model
    end
    V1=Array{Float64}[]
    for (k,v) in enumerate(V)
      append!(v,[1.0])
      push!(V1,vec((v')*transpose(affineMatrix)))
      pop!(V[k])
      pop!(V1[k])
    end
    if length(model)==2
      return V1,CV
    elseif length(model)==3
      return V1,CV,FV
    end
   end 
   return larApply0
end


#____________________plarApply______________________________

@everywhere function plarApply(affineMatrix)
  function plarApply0(model)
    if length(model)==2
      V,CV=deepcopy(model)
    elseif length(model)==3
      V,CV,FV = deepcopy(model)
    end
    V1=Array{Float64}[]
    V1=@sync @parallel (append!)for v in V
      append!(v,[1.0])
      [collect(vec((v')*transpose(affineMatrix)))]
    end
    for v in V1
      pop!(v)
    end
    if length(model)==2
      return fetch(V1),CV
    elseif length(model)==3
      return V1,CV,FV
    end
  end
  return plarApply0
end


#_________________________"larApply Tests"__________________________________

square=([[0,0],[0,1],[1,0],[1,1]],[[0,1,2,3]])
cubes=([[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,1,0],[1,1,1]],[[0,1,2,3,4,5,6,7]])

@testset "larApply Tests" begin
  @testset "2D" begin
    @testset "larApply Translation 2D" begin
      @test typeof(larApply(t(-0.5,-0.5))(square))==Tuple{Array{Array{Float64,N} where N,1},
      Array{Array{Int64,1},1}}
      @test larApply(t(-0.5,-0.5))(square)==([[-0.5,-0.5],[-0.5,0.5],[0.5,-0.5],[0.5,0.5]],
      [[0,1,2,3]])
    end
    
    @testset "larApply Scaling 2D" begin
      @test typeof(larApply(s(-0.5,-0.5))(square))==Tuple{Array{Array{Float64,N} where N,1},
      Array{Array{Int64,1},1}}
      @test larApply(s(-0.5,-0.5))(square)==([[0.0,0.0],[0.0,-0.5],[-0.5,0.0],[-0.5,-0.5]],
      [[0,1,2,3]])
    end
    
    @testset "larApply Rotation 2D" begin
      @test typeof(larApply(r(0))(square))==Tuple{Array{Array{Float64,N}  where N,1},
      Array{Array{Int64,1},1}}
      @test larApply(r(0))(square)==square
    end
  end
  
  @testset "3D" begin
    @testset "larApply Translation 3D" begin
      @test typeof(larApply(t(-0.5,-0.5,-0.5))(cubes))==Tuple{Array{Array{Float64,N}
      where N,1},Array{Array{Int64,1},1}}
      @test larApply(t(-0.5,-0.5,-0.5))(cubes)==([[-0.5,-0.5,-0.5],[-0.5,-0.5,0.5],
      [-0.5,0.5,-0.5],
      [-0.5,0.5,0.5],[0.5,-0.5,-0.5],[0.5,-0.5,0.5],[0.5,0.5,-0.5],[0.5,0.5,0.5]],
      [[0,1,2,3,4,5,6,7]])
    end
    
    @testset "larApply Scaling 3D" begin
      @test typeof(larApply(s(-0.5,-0.5,-0.5))(cubes))==Tuple{Array{Array{Float64,N} 
      where N,1},Array{Array{Int64,1},1}}
      @test larApply(s(-0.5,-0.5,-0.5))(cubes)==([[0.0,0.0,0.0],[0.0,0.0,-0.5],[0.0,-0.5,0.0],
      [0.0,-0.5,-0.5],[-0.5,0.0,0.0],[-0.5,0.0,-0.5],[-0.5,-0.5,0.0],[-0.5,-0.5,-0.5]],
      [[0,1,2,3,4,5,6,7]])
    end
    
    @testset "larApply Rotation 3D" begin
      @test typeof(larApply(r(pi,0,0))(cubes))==Tuple{Array{Array{Float64,N} where N,1},
      Array{Array{Int64,1},1}}
      @test isapprox(larApply(r(pi,0,0))(cubes)[1],[[0.0,0.0,0.0],[0.0,-1.22465e-16,-1.0],
      [0.0,-1.0,1.22465e-16],[0.0,-1.0,-1.0],[1.0,0.0,0.0],[1.0,-1.22465e-16,-1.0],
      [1.0,-1.0,1.22465e-16],[1.0,-1.0,-1.0]])
    end
  end
end

#_________________________"plarApply Tests"________________________________________

square=([[0,0],[0,1],[1,0],[1,1]],[[0,1,2,3]])
cubes=([[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,1,0],[1,1,1]],[[0,1,2,3,4,5,6,7]])
@testset "plarApply Tests" begin
  @testset "2D" begin
    @testset "plarApply Translation 2D" begin
      @test typeof(plarApply(t(-0.5,-0.5))(square))==Tuple{Array{Array{Float64,1},1},
      Array{Array{Int64,1},1}}
      @test plarApply(t(-0.5,-0.5))(square)==([[-0.5,-0.5],[-0.5,0.5],[0.5,-0.5],[0.5,0.5]],
      [[0,1,2,3]])
    end
    @testset "plarApply Scaling 2D" begin
      @test typeof(plarApply(s(-0.5,-0.5))(square))==Tuple{Array{Array{Float64,1},1},
      Array{Array{Int64,1},1}}
      @test plarApply(s(-0.5,-0.5))(square)==([[0.0,0.0],[0.0,-0.5],[-0.5,0.0],
      [-0.5,-0.5]],[[0,1,2,3]])
    end
    @testset "plarApply Rotation 2D" begin
      @test typeof(plarApply(r(0))(square))==Tuple{Array{Array{Float64,1},1},
      Array{Array{Int64,1},1}}
      @test plarApply(r(0))(square)==square
    end
  end
  @testset "3D" begin
    @testset "plarApply Translation 3D" begin
    @test typeof(plarApply(t(-0.5,-0.5,-0.5))(cubes))==Tuple{Array{Array{Float64,1},1},
    Array{Array{Int64,1},1}}
    @test plarApply(t(-0.5,-0.5,-0.5))(cubes)==([[-0.5,-0.5,-0.5],[-0.5,-0.5,0.5],[-0.5,0.5,-0.5],
    [-0.5,0.5,0.5],[0.5,-0.5,-0.5],[0.5,-0.5,0.5],[0.5,0.5,-0.5],[0.5,0.5,0.5]],[[0,1,2,3,4,5,6,7]])
  end 
    @testset "plarApply Scaling 3D" begin
      @test typeof(plarApply(s(-0.5,-0.5,-0.5))(cubes))==Tuple{Array{Array{Float64,1},1},
      Array{Array{Int64,1},1}}
      @test plarApply(s(-0.5,-0.5,-0.5))(cubes)==([[0.0,0.0,0.0],[0.0,0.0,-0.5],[0.0,-0.5,0.0],
      [0.0,-0.5,-0.5],[-0.5,0.0,0.0],[-0.5,0.0,-0.5],[-0.5,-0.5,0.0],[-0.5,-0.5,-0.5]],
      [[0,1,2,3,4,5,6,7]])
    end
    @testset "plarApply Rotation 3D" begin
      @test typeof(plarApply(r(pi,0.0,0.0))(cubes))==Tuple{Array{Array{Float64,1},1},
      Array{Array{Int64,1},1}}
      @test isapprox(plarApply(r(pi,0,0))(cubes)[1],[[0.0,0.0,0.0],[0.0,-1.22465e-16,-1.0],
      [0.0,-1.0,-1.22465e-16],[0.0,-1.0,-1.0],
      [1.0,0.0,0.0],[1.0,-1.22465e-16,-1.0],[1.0,-1.0,1.22465e-16],[1.0,-1.0,-1.0]])
    end
  end
end

#___________________________Execution time on PC_________________________________

times=[]
ptimes=[]
append!(ptimes,Time(plarApply(t(-0.5,-0.5)),[l[i]]) for i in range(1,length(l)))
append!(ptimes,Time(plarApply(t(-0.5,-0.5)),[l[i]]) for i in range(1,length(l)))

plot(times,xlabel="input",xlims=(0,length(times)+2),ylabel="time(s)",label=["Serial"])
plot(ptimes,xlabel="input",xlims=(0,length(times)+2),ylabel="time(s)",label=["Parallel"])

plot([times,ptimes],xlabel="input",xlims=(0,length(times)+2),ylabel="time(s)",
label=["Serial","Parallel"])
#____________________box__________________________
 
function box(model)
  if isa(model,Matrix)
    return []
  elseif isa(model,Struct)
    dummyModel=deepcopy(model)
    dummyModel.body=Any[]
    for term in model.body
      if isa(term,Struct)
	push!(dummyModel.body,[term.box,[0,1]])
      else
	push!(dummyModel.body,term)
      end
    end
    listOfModels=evalStruct(dummyModel)
    theMin,theMax=box(listOfModels[1])
    for theModel in listOfModels[2:end]
      modelMin,modelMax= box(theModel)
      for (k,val) in enumerate(modelMin)
	if val < theMin[k]
	  theMin[k]=val
	end
      end
      for (k,val) in enumerate(modelMax)
	if val > theMax[k]
	  theMax[k]=val
	end
      end
    end
    return Array[theMin,theMax]
  elseif (isa(model,Tuple)||isa(model,Array)) &&(length(model)==2||length(model)==3)
    V=model[1]
    theMin=[]
    theMax=[]
    for j in range(1,length(V[1]))
      Min=V[1][j]
      Max=V[1][j]
      for i in range(1,length(V))
	Min=min(Min,V[i][j])
	Max=max(Max,V[i][j])
      end
      push!(theMin,Min)
      push!(theMax,Max)
    end
    return Array[theMin,theMax]
  end
end

#____________________pbox____________________________

function pbox(model)
  if isa(model,Matrix)
    return []
  elseif isa(model,pStruct)
    dummyModel=deepcopy(model)
    dummyModel.body=Any[]
    @sync for term in model.body
      if isa(term,pStruct)
	push!(dummyModel.body,[term.box,[0,1]])
      else
	push!(dummyModel.body,term)
      end
    end
    listOfModels=pevalStruct(dummyModel)
    theMin,theMax=pbox(listOfModels[1])
    @sync for theModel in listOfModels[2:end]
      modelMin,modelMax= pbox(theModel)
      @async begin
	for (k,val) in enumerate(modelMin)
	  if (val < theMin[k])
	    theMin[k]=val
	  end
	end
	for (k,val) in enumerate(modelMax)
	  if (val > theMax[k])
	    theMax[k]=val
	  end
	end
      end
      end
      return Array[theMin,theMax]
  elseif (isa(model,Tuple) ||isa(model,Array))&&(length(model)==2 || length(model)==3)
    V=model[1]
    theMin=[]
    theMax=[]
    @sync for j in range(1,length(V[1]))
      Min=V[1][j]
      Max=V[1][j]
      for i in range(1,length(V))
	Min=min(Min,V[i][j])
	Max=max(Max,V[i][j])
      end
      @async begin
	push!(theMin,Min)
	push!(theMax,Max)
      end
    end
    return Array[theMin,theMax]
  end
end



#_________________________traversal_______________________________

function traversal(CTM,stack,obj,scene=[])
  for i in range(1,len(obj))
    if isa(obj.body[i],Matrix)
      CTM = CTM * obj.body[i]
    elseif (isa(obj.body[i],Tuple)||isa(obj.body[i],Array))&&
    (length(obj.body[i])==2||length(obj.body[i])==3)
      l=larApply(CTM)(obj.body[i])
      push!(scene,l)
    elseif isa(obj.body[i],Struct)
      push!(stack,CTM)
      traversal(CTM,stack,obj.body[i],scene)
      CTM = pop!(stack)
    end
  end
  return scene
end

#___________________ptraversal__________________________

@everywhere function ptraversal(CTM,stack,obj,scene=[])
  @sync for i in range(1,len(obj))
    if isa(obj.body[i],Matrix)
      CTM = CTM*obj.body[i]
    elseif (isa(obj.body[i],Tuple)||isa(obj.body[i],Array))&&
    (length(obj.body[i])==2||length(obj.body[i])==3)
      l=plarApply(CTM)(obj.body[i])
      push!(scene,l)
    elseif isa(obj.body[i],pStruct)
      push!(stack,CTM)
      ptraversal(CTM,stack,obj.body[i],scene)
      CTM = pop!(stack)
    end
  end
  return scene
end


#________________Struct________________________

#evalstruct

function evalStruct(self)
  dim = checkStruct(self.body)
  CTM, stack = eye(dim+1), []
  scene = traversal(CTM, stack, self, []) 
  return scene
end

#Struct:
  
type Struct
  body::Array
  box
  name::AbstractString
  dim
  category::AbstractString
  function Struct()
    self=new([],Nullable{Any},"new",Nullable{Any},"feature")
    self.name=string(object_id(self))
    return self
  end
  function Struct(data::Array)
    self=Struct()
    self.body=data
    self.box=box(self)
    self.dim=length(self.box[1])
    return self
  end
  function Struct(data::Array,name)
    self=Struct()
    self.body=[item for item in data]
    self.box=box(self)
    self.dim=length(self.box[1])
    self.name=string(name)
    return self
  end
  function Struct(data::Array,name,category)
    self=Struct()
    self.body=[item for item in data]
    self.box=box(self)
    self.dim=length(self.box[1])
    self.name=string(name)
    self.category=string(category)
    return self
  end
end
function name(self::Struct)
  return self.name
end
function category(self::Struct)
  return self.category
end
function len(self::Struct)
  return length(self.body)
end
function getitem(self::Struct,i::Int)
  return self.body[i]
end
function setitem(self::Struct,i,value)
  self.body[i]=value
end
function pprint(self::Struct)
  return "<Struct name: $(self.__name__())"
end
function set_name(self::Struct,name)
  self.name=string(name)
end
function clone(self::Struct,i=0)
    newObj=deepcopy(self)
    if i!=0
      newObj.name="$(self.__name__())_$(string(i))"
    end
    return newObj
end
function set_category(self::Struct,category)
  self.category=string(category)
end

#___________________pStruct______________________

#pevalStruct

function pevalStruct(self)
  dim = pcheckStruct(self.body)
  CTM, stack = eye(dim+1), []
  scene = ptraversal(CTM, stack, self, []) 
return scene
end

#pStruct

@everywhere type pStruct
  body::Array
  box
  name::AbstractString
  dim
  category::AbstractString
  function pStruct()
    self=new([],Nullable{Any},"new",Nullable{Any},"feature")
    self.name=string(object_id(self))
    return self
  end
  function pStruct(data::Array)
    self=pStruct()
    self.body=data
    self.box=pbox(self)
    self.dim=length(self.box[1])
    return self
  end
  function pStruct(data::Array,name)
    self=pStruct()
    self.body=[item for item in data]
    self.box=pbox(self)
    self.dim=length(self.box[1])
    self.name=string(name)
    return self
  end
  function pStruct(data::Array,name,category)
    self=pStruct()
    self.body=[item for item in data]
    self.box=pbox(self)
    self.dim=length(self.box[1])
    self.name=string(name)
    self.category=string(category)
    return self
  end
end
function name(self::pStruct)
  return self.name
end
function category(self::pStruct)
  return self.category
end
function len(self::pStruct)
  return length(self.body)
end
function getitem(self::pStruct,i::Int)
  return self.body[i]
end
function setitem(self::pStruct,i,value)
  self.body[i]=value
end
function pprint(self::pStruct)
  return "<Struct name: $(self.__name__())"
end
function set_name(self::pStruct,name)
  self.name=string(name)
end
function clone(self::pStruct,i=0)
  newObj=deepcopy(self)
  if i!=0
    newObj.name="$(self.__name__())_$(string(i))"
  end
  return newObj
end
function set_category(self::pStruct,category)
  self.category=string(category)
end

#___________________"Struct Tests"_________________________

square=([[0,0],[0,1],[1,0],[1,1]],[[0,1,2,3]])

@testset "Struct Tests" begin
	@test Struct([square]).body==[square]
	@test Struct([square]).dim==length(square[1][1])
	@test Struct([square]).box==[[0,0],[1,1]]
end

#__________________________"pStruct Tests"________________________

square=([[0,0],[0,1],[1,0],[1,1]],[[0,1,2,3]])

@testset "pStruct Tests" begin
    @test pStruct([square]).body==[square]
    @test pStruct([square]).dim==length(square[1][1])
    @test pStruct([square]).box==[[0,0],[1,1]]
    @test pStruct([square]).category=="feature"
    @test pStruct([square],"quadrato").name=="quadrato"
end


#______________________Execution time on PC__________________________________
times=[]
ptimes=[]

append!(times,Time(Struct,[[l[i]]]) for i in range(1,length(l)) )
append!(ptimes,Time(pStruct,[[l[i]]]) for i in range(1,length(l)-1))

plot(times,xlabel="input",xlims=(0,length(times)+2),ylabel="time(s)",label=["Serial"])
plot(ptimes,xlabel="input",xlims=(0,length(times)+2),ylabel="time(s)",label=["Parallel"])

plot([times,ptimes],xlabel="input",xlims=(0,length(times)+2),ylabel="time(s)",
label=["Serial","Parallel"])

#___________________Execution time on Tesla______________________________________

input=[1,10,20,50,10^2,2*10^2,5*10^2,10^3,2*10^3]
function timeStruct(model,input)
  t=Array{Float64}(length(input))
  pt=Array{Float64}(length(input))
  for i in range(1,length(input))
    structo=addn2D(input[i],model)
    Struct(structo.body)
    pStruct(structo.body)
    t[i]=@elapsed  Struct(structo.body)
    pt[i]=@elapsed pStruct(structo.body)
  end
  return t,pt
end

y,yp=timeStruct(square,input)
p=plot(y,xaxis="input",yaxis="time",xlims=(0,length(input)+1),
       ylims=(0,maximum(y)+0.5), label=["Serial"],lw=2)
       
pp=plot(yp,xaxis="input",yaxis="time",xlims=(0,length(input)+1),
       ylims=(0,maximum(y)+0.5),label=["Parallel"],lw=2)

yc=[y,yp]
pc=plot(yc,label=["Serial" "Parallel"])
###################################################################################################################################

#___________________"box Tests"_______________________

square=([[0,0],[0,1],[1,0],[1,1]],[[0,1,2,3]])

cubes=([[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,1,0],[1,1,1]],[[[0],[1],[2],[3],[4],[5],
[6],[7]],[[0,1],[2,3],[4,5],[6,7],[0,2],[1,3],[4,6],[5,7],[0,4],[1,5],[2,6],[3,7]],[[0,1,2,3],
[4,5,6,7],[0,1,4,5],[2,3,6,7],[0,2,4,6],[1,3,5,7]],[[0,1,2,3,4,5,6,7]]])

@testset "box Tests" begin
  @testset "box Tests 2D" begin
    @test typeof(box(square))==Array{Array,1}
    @test length(box(square))==2
    @test length(box(square)[1])==2
  end
  
  @testset "box Tests 3D" begin
    @test typeof(box(cubes))==Array{Array,1}
    @test length(box(cubes))==2
    @test length(box(cubes)[1])==3
  end
end

#_____________________"pbox Tests"___________________________________

square=([[0,0],[0,1],[1,0],[1,1]],[[0,1,2,3]])

cubes=([[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,1,0],[1,1,1]],[[[0],[1],[2],[3],[4],[5],
[6],[7]],[[0,1],[2,3],[4,5],[6,7],[0,2],[1,3],[4,6],[5,7],[0,4],[1,5],[2,6],[3,7]],[[0,1,2,3],[4,5,6,7],
[0,1,4,5],[2,3,6,7],[0,2,4,6],[1,3,5,7]],[[0,1,2,3,4,5,6,7]]])

@testset "pbox Tests" begin
  @testset "pbox Tests 2D" begin
    @test typeof(pbox(square))==Array{Array,1}
    @test length(pbox(square))==2
    @test length(pbox(square)[1])==2
  end
  
  @testset "pbox Tests 3D" begin
    @test typeof(pbox(cubes))==Array{Array,1}
    @test length(pbox(cubes))==2
    @test length(pbox(cubes)[1])==3
  end
end


#_______________________Execution time on PC_________________________

times=[]
ptimes=[]
input=[]

for i in range(1,length(l))
  push!(input,Struct([repeat([l[i]],outer=i)...]))
end
append!(times,Time(box,[input[i]]) for i in range(1,length(input)))
append!(ptimes,Time(pbox,[input[i]])  for i in range(1,length(input)))

plot(times,xlabel="input",xlims=(0,length(times)+2),ylabel="time(s)",label=["Serial"])
plot(ptimes,xlabel="input",xlims=(0,length(times)+2),ylabel="time(s)",label=["Parallel"])

plot([times,ptimes],xlabel="input",xlims=(0,length(ptimes)+2),ylabel="time(s)",
label=["Serial","Parallel"])
#_____________________Execution time on Tesla_______________________________

input=[1,10,20,50,10^2,2*10^2,5*10^2,10^3,2*10^3]
function timeFstruct(f::Function,pf::Function,model,input)
  t=Array{Float64}(length(input))
  pt=Array{Float64}(length(input))
  for i in range(1,length(input))
    structo=addn2D(input[i],model)
    pstructo=pStruct(structo.body)
    f(structo)
    pf(pstructo)
    t[i]=@elapsed f(structo)
    pt[i]=@elapsed pf(pstructo)
  end
  return t,pt
end
y,yp=timeFstruct(box,pbox,square,input)
p=plot(y,xaxis="input",yaxis="time",xlims=(0,length(input)+1),ylims=(0,maximum(y)+0.5),
label=["Serial"],lw=2)
pp=plot(yp,xaxis="input",yaxis="time",xlims=(0,length(input)+1),ylims=(0,maximum(y)+0.5),
label=["Parallel"],lw=2)

yc=[y,yp]
pc=plot(y,xaxis="input",yaxis="time",xlims=(0,length(input)+1),ylims=(0,maximum(y)+0.5),
label=["Serial" "Parallel"],lw=2)

####################################################################################################################################

#_______________________"traversal Tests"_______________________

@testset "traversal" begin
   square=([[0, 0], [0, 1], [1, 0], [1, 1]], [[0, 1, 2, 3]])
   
   @everywhere structure=Struct([square])
   @everywhere d=checkStruct(structure.body)
   @test length(traversal(eye(d+1),[],structure,[]))==length(structure.body)
   @test typeof(traversal(eye(d+1),[],structure,[]))==Array{Any,1}
end

#______________________"ptraversal tests"____________________________

square=([[0, 0], [0, 1], [1, 0], [1, 1]], [[0, 1, 2, 3]])
structure=pStruct([square])
d=pcheckStruct(structure.body)

@testset "ptraversal" begin
  @test length(ptraversal(eye(d+1),[],structure,[]))==length(structure.body)
  @test typeof(ptraversal(eye(d+1),[],structure,[]))==Array{Any,1}
end

#__________________Execution time on PC______________________

times=[]
ptimes=[]
input=[]
dim=[]

for i in range(1,length(l)-1)
  push!(input,Struct([repeat([l[i]],outer=i)...]))
end

append!(dim,checkStruct(input[i].body) for i in range(1,length(input)))

for i in range(1,length(input))
  args=(eye(dim[i]+1),[],input[i],[])
  append!(times,Time(traversal,(args...)))
end
for i in range(1,length(input))
  args=(eye(dim[i]+1),[],input[i],[])
  append!(ptimes,Time(ptraversal,(args...)))
end

plot(times,xlabel="input",xlims=(0,length(times)+2),ylabel="time(s)",label=["Serial"])
plot(ptimes,xlabel="input",xlims=(0,length(times)+2),ylabel="time(s)",label=["Parallel"])


plot([times,ptimes],xlabel="input",xlims=(0,length(ptimes)+2),ylabel="time(s)",
label=["Serial","Parallel"])

#_____________________Execution time on Tesla____________________________

input=[1,10,20,50,10^2,2*10^2,5*10^2,10^3,2*10^3]

function timeTraversal(model,input)
   t=Array{Float64}(length(input))
   pt=Array{Float64}(length(input))
   for i in range(1,length(input))
      structo=addn2D(input[i],model)
      pstructo=pStruct(structo.body)
      dim=checkStruct(structo.body)
      traversal(eye(dim+1),[],structo,[])
      ptraversal(eye(dim+1),[],pstructo,[])
      pt[i]=@elapsed ptraversal(eye(dim+1),[],pstructo,[])
      t[i]=@elapsed traversal(eye(dim+1),[],structo,[])
   end
   return t,pt
end

y,yp=timeTraversal(square,input)


p=plot(y,xaxis="input",yaxis="time",xlims=(0,length(input)+1),ylims=(0,maximum(y)+0.5),
  label=["Serial"],lw=2)  
pp=plot(yp,xaxis="input",yaxis="time",xlims=(0,length(input)+1),
       ylims=(0,maximum(y)+0.5),label=["Parallel"],lw=2)

yc=[y,yp]
pc=plot(yc,label=["Serial" "Parallel"])


####################################################################################################################################

#______________________________embedTraversal_______________________

  
function embedTraversal(cloned,obj,n,suffix)
  for i in range(1,len(obj))
    if isa(obj.body[i],Matrix)
      mat=obj.body[i]
      d,d=size(mat)
      newMat=eye(d+n*1)
      for h in range(1,d-1)
	for k in range(1,d-1)
	  newMat[h,k]=mat[h,k]
	end
	newMat[h,d-1+n*1]=mat[h,d-1]
      end
      append!(cloned.body,newMat)
    elseif (isa(obj.body[i],Tuple) ||isa(obj.body[i],Array))&&length(obj.body[i])==3 
      V,FV,EV=obj.body[i]
      dimadd=fill([0.0],n)
      for k in dimadd
	for v in V
	  append!(v,k)
	end
      end
      append!(cloned.body,[(V,FV,EV)])
    elseif (isa(obj.body[i],Tuple) ||isa(obj.body[i],Array))&&length(obj.body[i])==2
      V,EV=deepcopy(obj.body[i])
      dimadd=fill([0.0],n)
      for k in dimadd
	for v in V
	  append!(v,k)
	end
      end
      append!(cloned.body,[(V,EV)])
    elseif isa(obj.body[i],Struct)
      newObj=Struct()
      newObj.box=hcat((obj.body[i].box,[fill([0],n),fill([0],n)]))
      newObj.category=obj.body[i].category
      append!(cloned.body,embedTraversal(newObj,obj.body[i],n,suffix))
    end
  end
  return cloned
end

#__________________________pembedTraversal_______________________________

@everywhere function pembedTraversal(cloned,obj,n,suffix)
  for i in range(1,len(obj))
    if (isa(obj.body[i],Matrix) || isa(obj.body[i],SharedArray))
      mat=obj.body[i]
      d,d=size(mat)
      newMat=eye(d+n*1)
      @sync begin
	for h in range(1,d-1)
	  @async begin
	    for k in range(1,d-1)
	      newMat[h,k]=mat[h,k]
	    end
	  end
	  newMat[h,d-1+n*1]=mat[h,d-1]
	end
      end
      append!(cloned.body,newMat)
    elseif (isa(obj.body[i],Tuple) ||isa(obj.body[i],Array))&& length(obj.body[i])==2 
      V,EV=deepcopy(obj.body[i])
      dimadd=fill([0.0],n)
      @sync begin
	for k in dimadd
	  @async begin
	    for v in V
	      append!(v,k)
	    end
	  end
	end
      end
      append!(cloned.body,[(V,EV)])
    elseif (isa(obj.body[i],Tuple) ||isa(obj.body[i],Array))&&length(obj.body[i])==3
      V,FV,EV=obj.body[i]
      dimadd=fill([0.0],n)
      @sync begin
	for k in dimadd
	  @async begin
	    for v in V
	      append!(v,k)
	    end
	  end
	end
      end
      append!(cloned.body,[(V,FV,EV)])
    elseif isa(obj.body[i],pStruct)
      newObj=pStruct()
      @async begin
	newObj.box=hcat((obj.body[i].box,[fill([0],n),fill([0],n)]))
	newObj.category=obj.body[i].category
	append!(cloned.body,embedTraversal(newObj,obj.body[i],n,suffix))
      end
    end
  end
  return cloned
end

#__________________"embedTraversal Tests"_________________________________

square=([[0,0],[0,1],[1,0],[1,1]],[[0,1,2,3]])
x=Struct([square])

@testset "embedTraversal Tests" begin
  @test length(embedTraversal(deepcopy(x),deepcopy(x),1,"New").body[2][1][1])==
  length(x.body[1][1][1])+1
  #in this case n=1, but generally:
  # length(length(embedTraversal(x,x,1,"New")=length(x.body[1][1][1])+n
  @test length(embedTraversal(deepcopy(x),deepcopy(x),3,"New").body[2][1][1])==
  length(x.body[1][1][1])+3
  @test typeof(embedTraversal(deepcopy(x),deepcopy(x),1,"New"))==Struct	
end


#____________________________"pembedTraversal Tests"__________________________

square=([[0,0],[0,1],[1,0],[1,1]],[[0,1,2,3]])
x=pStruct([square])

@testset "pembedTraversal Tests" begin
  @test length(pembedTraversal(deepcopy(x),deepcopy(x),1,"New").body[2][1][1])==
  length(x.body[1][1][1])+1
  #in this case n=1, but generally:
  #length(length(embedTraversal(x,x,1,"New")=length(x.body[1][1][1])+n
  @test length(pembedTraversal(deepcopy(x),deepcopy(x),3,"New").body[2][1][1])==
  length(x.body[1][1][1])+3
  @test typeof(pembedTraversal(deepcopy(x),deepcopy(x),1,"New"))==pStruct	
end

#_________________________Execution time on Tesla________________________________________

input=[1,10,20,50,10^2,2*10^2,5*10^2,10^3,2*10^3]
function timeEmbedTraversal(model,input)
  t=Array{Float64}(length(input))
  pt=Array{Float64}(length(input))
  for i in range(1,length(input))
    structo=addn2D(input[i],model)
    pstructo=pStruct(structo.body)
    cloned=Struct()
    cloned.box=hcat((structo.box,[fill([0],10),fill([0],10)]))
    cloned.name=string(object_id(cloned))
    cloned.category=structo.category
    cloned.dim=structo.dim+10
    pcloned=pStruct()
    pcloned.box=hcat((pstructo.box,[fill([0],10),fill([0],10)]))
    pcloned.name=string(object_id(pcloned))
    pcloned.category=pstructo.category
    pcloned.dim=pstructo.dim+10
    embedTraversal(cloned,structo,10,"New")
    pembedTraversal(pcloned,pstructo,10,"New")
    t[i]=@elapsed embedTraversal(cloned,structo,10,"New")
    pt[i]=@elapsed pembedTraversal(pcloned,pstructo,10,"New")
  end
  return t,pt
end
y,yp=timeEmbedTraversal(square,input)
p=plot(input,y,xaxis="input",yaxis="time",xlims=(0,length(input)+1),
       ylims=(0,maximum(y)+0.5), label=["Serial"],lw=2)       
pp=plot(input,yp,xaxis="input",yaxis="time",xlims=(0,length(input)+1),
       ylims=(0,maximum(y)+0.5),label=["Parallel"],lw=2)      


yc=[y,yp]
pc=plot(input,yc,label=["Serial" "Parallel"])

#__________________________embedStruct____________________________________

function embedStruct(n)
  function embedStruct0(self,suffix="New")
    if n==0
      return self, length(self.box[1])
    end
    cloned=Struct()
    cloned.box=hcat((self.box,[fill([0],n),fill([0],n)]))
    cloned.name=string(object_id(cloned))
    cloned.category=self.category
    cloned.dim=self.dim+n
    cloned=embedTraversal(cloned,self,n,suffix)
    return cloned
 end
 return embedStruct0
end

#___________________________pembedStruct__________________________

@everywhere function pembedStruct(n)
  function pembedStruct0(self,suffix="New")
    if n==0
      return self, length(self.box[1])
    end
    cloned=pStruct()
    cloned.box=hcat((self.box,[fill([0],n),fill([0],n)]))
    cloned.name=string(object_id(cloned))
    cloned.category=self.category
    cloned.dim=self.dim+n
    cloned=pembedTraversal(cloned,self,n,suffix)
    return cloned
  end
  return pembedStruct0
end

#______________________"embedStruct Tests"____________________________

square=([[0,0],[0,1],[1,0],[1,1]],[[0,1,2,3]])
x=Struct([square])	

@testset "embedStruct Tests" begin
  @test length(embedStruct(1)(x).body[1][1][1])==length(x.body[1][1][1])+1
  #in this case n = 1, but generally: 
  #length(embedStruct(n)(x).body[1][1][1])=length(x.body[1][1][1])+n
  @test length(embedStruct(3)(x).body[1][1][1])==length(x.body[1][1][1])+3
  @test typeof(embedStruct(1)(x))==Struct	
end

#__________________"pembedStruct Tests"________________________________

square=([[0,0],[0,1],[1,0],[1,1]],[[0,1,2,3]])
x=pStruct([square])

@testset "pembedStruct Tests" begin
  @test length(pembedStruct(1)(x).body[1][1][1])==length(x.body[1][1][1])+1 
  #in this case n = 1, but generally: 
  #length(embedStruct(n)(x).body[1][1][1])=length(x.body[1][1][1])+n
  @test length(pembedStruct(3)(x).body[1][1][1])==length(x.body[1][1][1])+3
  @test typeof(pembedStruct(1)(x))==pStruct	
end

#______________________________Execution time on Tesla____________________________
input=[1,10,20,50,10^2,2*10^2,5*10^2,10^3,2*10^3]

function timeEmbedStruct(n,model,input)
  t=Array{Float64}(length(input))
  pt=Array{Float64}(length(input))
  for i in range(1,length(input))
    structo=addn2D(input[i],model)
    pstructo=pStruct(structo.body)
    embedStruct(n)(structo)
    pembedStruct(n)(pstructo)
    t[i]=@elapsed embedStruct(n)(structo)
    pt[i]=@elapsed pembedStruct(n)(pstructo)
  end
  return t,pt
end
y,yp=timeEmbedStruct(10,square,input)

p=plot(input,y,xaxis="input",yaxis="time",xlims=(0,length(input)+1),
       ylims=(0,maximum(y)+0.5), label=["Serial"],lw=2)       
pp=plot(input,yp,xaxis="input",yaxis="time",xlims=(0,length(input)+1),
       ylims=(0,maximum(y)+0.5),label=["Parallel"],lw=2)

yc=[y,yp]
pc=plot(input,yc,label=["Serial" "Parallel"])

#______________removeDups_________________________________

function removeDups(CW)
  CW=collect(Set(CW))
  CWs=collect(map(sort,CW))
  no_duplicates=Dict()
  for f in CWs
    no_duplicates[f] = []
  end
  for f in CW
    no_duplicates[sort(f)]=[f]
  end
  CW=[f[1] for f in values(no_duplicates)]
  return CW
end 


#______________________________premoveDups___________________________

function premoveDups(CW)
  CW=collect(Set(CW))
  CWs=collect(@sync pmap(sort,CW))
  no_duplicates=Dict()
  @parallel for f in CWs
    no_duplicates[f] = []
  end
  @parallel for f in CW
    no_duplicates[sort(f)]=[f]
  end
  @parallel for f in values(no_duplicates)
    append!(CW,f[1])
  end
  return CW
end

#______________________"removeDups Tests"____________________________

CW1=[[0,1,2,3],[4,5,6,7],[0,1,4,5],[2,3,6,7],[0,2,4,6],[1,3,5,7],[4,5,6,7],[8,9,10,11],
[4,5,8,9],[6,7,10,11],[4,6,8,10],[5,7,9,11]]
CW2=[[0,1,2,3],[4,5,6,7],[8,9,10,11],[12,13,14,15],[16,17,18,19]]

@testset "removeDups Tests" begin
  @testset "removeDups 3D" begin
    @test length(removeDups(CW1))<= length(CW1)
      @test typeof(removeDups(CW1))==Array{Array{Int64,1},1}
   end
   @testset "removeDups 2D" begin
      @test length(removeDups(CW2))<= length(CW2)
      @test typeof(removeDups(CW2))==Array{Array{Int64,1},1}
   end
end

#____________________________"premoveDups Tests"_______________________________________

CW1=[[0,1,2,3],[4,5,6,7],[0,1,4,5],[2,3,6,7],[0,2,4,6],[1,3,5,7],[4,5,6,7],[8,9,10,11],
[4,5,8,9],[6,7,10,11],[4,6,8,10],[5,7,9,11]]
CW2=[[0,1,2,3],[4,5,6,7],[8,9,10,11],[12,13,14,15],[16,17,18,19]]

@testset "premoveDups Tests" begin
   @testset "premoveDups 3D" begin
      @test length(premoveDups(CW1))<= length(CW1)
      @test typeof(premoveDups(CW1))==Array{Array{Int64,1},1}
   end
   @testset "premoveDups 2D" begin
      @test length(premoveDups(CW2))<= length(CW2)
      @test typeof(premoveDups(CW2))==Array{Array{Int64,1},1}
    end
end

#__________________________struct2lar______________________________

#fixedPrec

function fixedPrec(PRECISION)
  function fixedPrec0(value) 
    out=round.(value,PRECISION)
    if out==-0.0
      out=0.0
    end
    return string(out)
  end
  return fixedPrec0
end

#vcode 

function vcode(PRECISION=4)
  function vcode0(vect)
    return fixedPrec(PRECISION)(vect)
  end
  return vcode0
end

#struct2lar

function struct2lar(structure)
  listOfModels=evalStruct(structure)
  vertDict= Dict()
  index,defaultValue,CW,W,FW = -1,-1,[],[],[]
  for model in listOfModels
    if  length(model)==2
      V,FV=model
    elseif lenght(model)==3
      V,FV,EV=model
    end
    for (k,incell) in enumerate(FV)
      outcell=[ ]
      for v in incell
	key=vcode(4)(V[v+1])
	if get(vertDict,key,defaultValue)==defaultValue
	  index =index+1
	  vertDict[key]=index
	  append!(outcell,index)
	  append!(W,[eval(parse(key))])
	else
	  append!(outcell,vertDict[key])
	end
      end
      append!(CW,[outcell])
    end
    if length(model)==3
      for (k,incell) in enumerate(FV)
	outcell=[]
	for v in incell
	  key=vcode(4)(V[v+1])
	  if get(vertDict,key,defaultValue)==defaultValue
	    index =index+1
	    vertDict[key]=index
	    append!(outcell,[index])
	    append!(W,[eval(parse(key))])
	  else
	    append!(outcell,vertDict[key])
	  end
	end
	append!(FW,[outcell])
      end
    end
  end
  if length(listOfModels[end])==2
    if length(CW[1])==2
      CW=map(Tuple,map(sort,CW))
    else
      CW=removeDups(CW)
    end
    return W,CW
  end
  if length(listOfModels[end])==3
    FW=map(Tuple,map(sort,FW))
    CW=removeDups(CW)
    return W,CW,FW
    end
end

#___________________________pstruct2lar______________________________

@everywhere function pstruct2lar(structure)
  listOfModels=pevalStruct(structure)
  vertDict= Dict()
  index,defaultValue,CW,W,FW = -1,-1,[],[],[]
  for model in listOfModels
    if  length(model)==2
      V,FV=model
    elseif lenght(model)==3
      V,FV,EV=model
    end
    @sync begin
      for (k,incell) in enumerate(FV)
	outcell=[]
	@async begin
	  for v in incell
	    key=vcode(4)(V[v+1])
	    if get(vertDict,key,defaultValue)==defaultValue
	      index =index+1
	      vertDict[key]=index
	      append!(outcell,index)
	      append!(W,[eval(parse(key))])
	    else
	      append!(outcell,vertDict[key])
	    end
	  end
	end
	append!(CW,[outcell])
      end
    end
    if length(model)==3
      @sync begin
	for (k,incell) in enumerate(FV)
	  outcell=[]
	  @async begin
	    for v in incell
	      key=vcode(4)(V[v+1])
	      if get(vertDict,key,defaultValue)==defaultValue
		index =index+1
		vertDict[key]=index
		append!(outcell,[index])
		append!(W,[eval(parse(key))])
	      else
		append!(outcell,vertDict[key])
	      end
	    end
	  end
	  append!(FW,[outcell])
	end
      end
    end
  end
  if length(listOfModels[end])==2
    if length(CW[1])==2
      CW=pmap(Tuple,pmap(sort,CW))
    else
      CW=premoveDups(CW)
    end
    return W,CW
  end
  if length(listOfModels[end])==3
    FW=pmap(Tuple,pmap(sort,FW)) 
    CW=premoveDups(CW)
    return W,CW,FW
  end
end

#_______________________"struct2lar Tests"________________________

@testset "struct2lar" begin
  @testset "struct2lar 2D" begin
    square=([[0, 0],[0,1],[1,0],[1,1]],[[0,1,2,3]])
    table=larApply(t(-0.5,-0.5))(square)
    structure=Struct([repeat([table,r(pi/2)],outer=2)...])
    @test typeof(struct2lar(structure))==Tuple{Array{Any,1},Array{Array{Any,1},1}}
    @test length(struct2lar(structure)[1][1])==2
  end
  
  @testset "struct2lar 3D" begin
    BV=[[0,1,2,3],[4,5,6,7],[0,1,4,5],[2,3,6,7],[0,2,4,6],[1,3,5,7]]
    V=[[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,1,0],[1,1,1]]
    block=[V,BV]
    structure=Struct(repeat([block,t(1,0,0)],outer=2));
    @test typeof(struct2lar(structure))==Tuple{Array{Any,1},Array{Array{Any,1},1}}
    @test length(struct2lar(structure)[1][1])==3
  end
end

#___________________________"pstruct2lar Tests"__________________________________

@testset "pstruct2lar" begin
  @testset "pstruct2lar 2D" begin
    square=([[0,0],[0,1],[1,0],[1,1]],[[0,1,2,3]])
    table=plarApply(t(-0.5,-0.5))(square)
    structure=pStruct([repeat([table,r(pi/2)],outer=2)...])
    @test typeof(pstruct2lar(structure))==Tuple{Array{Any,1},Array{Any,1}}
    @test length(pstruct2lar(structure)[1][1])==2
  end
  
  @testset "pstruct2lar 3D" begin
    BV=[[0,1,2,3],[4,5,6,7],[0,1,4,5],[2,3,6,7],[0,2,4,6],[1,3,5,7]]
    V=[[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,1,0],[1,1,1]]
    block=[V,BV]
    structure=pStruct(repeat([block,t(1,0,0)],outer=2));
    @test typeof(pstruct2lar(structure))==Tuple{Array{Any,1},Array{Any,1}}
    @test length(pstruct2lar(structure)[1][1])==3
   end
end

#_______________Execution time on PC_____________________________
times=[]
ptimes=[]
input=[]

for i in range(1,length(l)-1)
  push!(input,Struct([repeat([l[i]],outer=i)...]))
end
append!(times,Time(struct2lar,[input[i]]) for i in range(1,length(input)))
append!(ptimes,Time(pstruct2lar,[input[i]]) for i in range(1,length(input)))


plot(times,xlabel="input",xlims=(0,length(times)+2),ylabel="time(s)",label=["Serial"])
plot(ptimes,xlabel="input",xlims=(0,length(times)+2),ylabel="time(s)",label=["Parallel"])

plot([times,ptimes],xlabel="input",xlims=(0,length(times)+2),ylabel="time(s)",
label=["Serial","Parallel"])

#________________________Execution time on Tesla_______________________________________

input=[1,10,20,50,10^2,2*10^2,5*10^2,10^3,2*10^3]
y,yp=timeFstruct(struct2lar,pstruct2lar,square,input)

p=plot(input,y,xaxis="input",yaxis="time",xlims=(0,length(input)+1),
     ylims=(0,maximum(y)+0.5), label=["Serial"],lw=2)
pp=plot(input,yp,xaxis="input",yaxis="time",xlims=(0,length(input)+1),
       ylims=(0,maximum(y)+0.5),label=["Parallel"],lw=2)

yc=[y,yp]
pc=plot(input,yc,label=["Serial" "Parallel"])

#_____________________larRemoveVertices_______________________

function larRemoveVertices(V,FV)
  vertDict= Dict()
  index,defaultValue,CW,W,FW = -1,-1,[],[],[]
  for (k,incell) in enumerate(FV)
    outcell=[]
    for v in incell
      key=vcode(4)(V[v+1])
      if get(vertDict,key,defaultValue)==defaultValue
	index =index+1
	vertDict[key]=index
	append!(outcell,index)
	append!(W,[eval(parse(key))])
      else
	append!(outcell,vertDict[key])
      end
    end
    append!(FW,[outcell])
  end
  return W,FW
end

#_____________________________plarRemoveVertices___________________________

@everywhere function plarRemoveVertices(V,FV)
  vertDict= Dict()
  index,defaultValue,CW,W,FW = -1,-1,[],[],[]
  @async begin
    for (k,incell) in enumerate(FV)
      outcell=[]
      @sync begin
	for v in incell
	  key=vcode(4)(V[v+1])
	  if get(vertDict,key,defaultValue)==defaultValue
	    index =index+1
	    vertDict[key]=index
	    append!(outcell,index)
	    append!(W,[eval(parse(key))])
	  else
	    append!(outcell,vertDict[key])
	  end
	end
      end
      append!(FW,[outcell])
    end
  end
  return W,FW
end

#_______________________"larRemoveVertices Tests"______________________________

V=[[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1, 1, 0], [1, 1, 1]]
FV=[[0,1,2,3],[4,5,6,7],[0,1,4,5],[2,3,6,7],[0,2,4,6],[1,3,5,7]]

@testset "larRemoveVertices Tests" begin
  @test typeof(larRemoveVertices(V,FV))==Tuple{Array{Any,1},Array{Any,1}}
  @test length(larRemoveVertices(V,FV)[1])<= length(V)
end

#_______________________"plarRemoveVertices Tests"_________________________________

V=[[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1, 1, 0], [1, 1, 1]]
FV=[[0,1,2,3],[4,5,6,7],[0,1,4,5],[2,3,6,7],[0,2,4,6],[1,3,5,7]]

@testset "plarRemoveVertices Tests" begin
  @test typeof(plarRemoveVertices(V,FV))==Tuple{Array{Any,1},Array{Any,1}}
  @test length(plarRemoveVertices(V,FV)[1])<= length(V)
end

#__________________Execution time on PC__________________________________
times=[]
ptimes=[]
append!(times,Time(larRemoveVertices,[l[i][1],l[i][2]]) for i in range(1,length(l)-1))
append!(ptimes,Time(plarRemoveVertices,[l[i][1],l[i][2]]) for i in range(1,length(l)-1))


plot(times,xlabel="input",xlims=(0,length(times)+2),ylabel="time(s)",label=["Serial"])
plot(ptimes,xlabel="input",xlims=(0,length(times)+2),ylabel="time(s)",label=["Parallel"])

plot([times,ptimes],xlabel="input",xlims=(0,length(times)+2),ylabel="time(s)",
label=["Serial","Parallel"])

#________________________Examples____________________________________

#Examples 1

square= ([[0,0],[0,1],[1,0],[1,1]],[[0,1,2,3]])
table= larApply(t(-.5,-.5))(square)
chair= larApply(s(.35,.35))(table)
chair1= larApply(t(.75,0))(chair)
chair2= larApply(r(pi/2))(chair1)
chair3= larApply(r(pi/2))(chair2)
chair4= larApply(r(pi/2))(chair3)

#Examples 2

desk=([[0, 0],[0, 1],[0, 2],[1, 0],[1, 1],[1, 2],[2, 0],[2, 1],[2, 2],[3, 0],[3, 1],[3, 2],
[4, 0],[4, 1],[4, 2],[5, 0],[5, 1],[5, 2],[6, 0],[6, 1],[6, 2],[7, 0],[7, 1],[7, 2],[8, 0],
[8, 1],[8, 2],[9, 0],[9, 1],[9, 2],[10, 0],[10, 1],[10, 2],[11, 0],[11, 1],[11, 2],[12, 0],
[12, 1],[12, 2],[13, 0],[13, 1],[13, 2],[14, 0],[14, 1],[14, 2],[15, 0],[15, 1],[15, 2],
[16, 0],[16, 1],[16, 2],[17, 0],[17, 1],[17, 2],[18, 0],[18, 1],[18, 2],[19, 0],[19, 1],
[19, 2],[20, 0],[20, 1],[20, 2]], [[0, 1, 3, 4],[1, 2, 4, 5],[3, 4, 6, 7],[4, 5, 7, 8],
[6, 7, 9, 10],[7, 8, 10, 11],[9, 10, 12, 13],[10, 11, 13, 14],[12, 13, 15, 16],
[13, 14, 16, 17],[15, 16, 18, 19],[16, 17, 19, 20],[18, 19, 21, 22],[19, 20, 22, 23],
[21, 22, 24, 25],[22, 23, 25, 26],[24, 25, 27, 28],[25, 26, 28, 29],[27, 28, 30, 31],
[28, 29, 31, 32],[30, 31, 33, 34],[31, 32, 34, 35],[33, 34, 36, 37],[34, 35, 37, 38],
[36, 37, 39, 40],[37, 38, 40, 41],[39, 40, 42, 43],[40, 41, 43, 44],[42, 43, 45, 46],
[43, 44, 46, 47],[45, 46, 48, 49],[46, 47, 49, 50],[48, 49, 51, 52],[49, 50, 52, 53],
[51, 52, 54, 55],[52, 53, 55, 56],[54, 55, 57, 58],[55, 56, 58, 59],[57, 58, 60, 61],
[58, 59, 61, 62]])
seat=larApply(s(0.02,0.3))(desk)
seat=larApply(t(0.55,-0.85))(seat)
line=Struct([desk,repeat([seat,t(1.65,0)],outer=12)...])
line=Struct([line,t(23,0),line])
lines=Struct([repeat([line,t(0,-3)],outer=6)...])
teacherdesk=larApply(s(0.5,0.8))(desk)
teacherdesk=larApply(t(16.5,5))(teacherdesk)
chair=larApply(s(0.15,0.5))(teacherdesk)
chair=larApply(t(18,5))(chair)
classroom=Struct([teacherdesk,chair,lines])
class=evalStruct(classroom)

