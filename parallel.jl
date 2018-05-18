function prepKey(args)
               v=join(args,",")
return(v)
end


#____________________________________________________________________________________________________________________________

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

#____________________________________________________________________________________________________________________________

function vcode(PRECISION=4)
	function vcode0(vect)
		#return prepKey(map(fixedPrec(PRECISION),vect)) dovrebbe essere così ma in julia non serve mettere prepKey per ottenere lo stesso risultato di python probabilmente perchè il map funziona in maniera differente
		return fixedPrec(PRECISION)(vect) #quando vect è contiene più array bisogna usare map
	end
	return vcode0
end

#____________________________________________________________________________________________________________________________


@everywhere function pt(args...)
	d=length(args)
	mat=SharedArray(eye(d+1))
	@parallel for i in range(1,d)
		mat[i,d+1]=args[i]
	end
	return mat
end

#____________________________________________________________________________________________________________________________

@everywhere function ps(args...)
	d=length(args)
	mat=SharedArray(eye(d+1))
	@parallel for k in range(1,d)
		mat[k,k]=args[k]
	end
	return mat
end

#____________________________________________________________________________________________________________________________


@everywhere function pr(args...)
    args = collect(args)
    n = length(args)
	mat=eye(3) #vedere se farla diventare una sharedArray
    if n == 1 # rotation in 2D
	@async begin
		angle = args[1]
		#COS=remotecall_fetch(cos,2,angle)
			#COS=remotecall(cos,2,angle)
		COS = cos(angle)
		SIN = sin(angle)
		@sync begin
			
        		mat[1,1] = COS
			mat[1,2] = -SIN
        		mat[2,1] = SIN
			mat[2,2] = COS
		end
	
   end
   end
	if n == 3 # rotation in 3D
        	@async begin
			angle = norm(args); axis = normalize(args)
        		COS = cos(angle); SIN= sin(angle)
			@sync begin
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
		end
	end
	return mat
end

#____________________________________________________________________________________________________________________________

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

#____________________________________________________________________________________________________________________________

 
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


#____________________________________________________________________________________________________________________________


#classe Struct



@everywhere type Struct
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
	function print(self::Struct)
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

#____________________________________________________________________________________________________________________________

@everywhere function plarBoundary(self)
	data=pstruct2lar(self)
	if length(data)==3
		V,FV,EV=data
		return V,FV,EV
	else
		return "<Struct name $(self.__name__()) : boundary non computable"
	end
end


#____________________________________________________________________________________________________________________________


@everywhere function pstruct2lar(structure)
	listOfModels=evalStruct(structure)
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
			CW=pmap(Tuple,map(sort,CW))
		else
			CW=premoveDups(CW)
		end
		return W,CW
	end
	
	if length(listOfModels[end])==3
		FW=map(Tuple,map(sort,FW))
		CW=premoveDups(CW)
		return W,CW,FW
	end
end



#____________________________________________________________________________________________________________________________


function embedTraversal(cloned,obj,n,suffix)
	for i in range(1,length(obj))
		if (isa(obj[i],tuple) ||isa(obj[i],Array))&& length(obj[i])==2 #potrebbe dare problemi in ndarray
			V,EV=obj[i]
			#V=[append!(v,fill(0.0,n)) for v in V]    #provare a vedere se girano entrambe allo stesso modo
			dimadd=fill([0.0],n)
			for k in dimadd
				for v in V
					append!(v,k)
				end
			end
			append!(cloned.body,(V,EV))
		elseif (isa(obj[i],tuple) ||isa(obj[i],Array))&& length(obj[i])==3 #potrebbe dare problemi in ndarray
			V,FV,EV=obj[i]
			dimadd=fill([0.0],n)
			for k in dimadd
				for v in V
					append!(v,k)
				end
			end
			append!(cloned.body,(V,FV,EV))
		elseif isa(obj[i],Matrix)
			mat=obj[i]
			d,d=size(mat)
			newMat=eye(d+n*1)
			for h in range(1,d-1)
				for k in range(1,d-1)
					newMat[h,k]=mat[h,k]
				end
				newMat[h,d-1+n*1]=mat[h,d-1]
			end
			append!(cloned.body,newMat) 
		elseif isa(obj[i],Struct)
			newObj=Struct()
			newObj.box=hcat((obj[i].box,[fill([0],n),fill([0],n)]))
			newObj.category=obj[i].category
			append!(cloned.body,embedTraversal(newObj,obj[i],n,suffix))
		end
	end
	return cloned
end

#____________________________________________________________________________________________________________________________


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

#____________________________________________________________________________________________________________________________

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
		#dim=checkStruct(listOfModels)
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

	elseif (isa(model,Tuple) ||isa(model,Array))&& (length(model)==2 || length(model)==3)
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

 

 

#____________________________________________________________________________________________________________________________


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

#____________________________________________________________________________________________________________________________

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



#____________________________________________________________________________________________________________________________


function traversal(CTM,stack,obj,scene=[])
	for i in range(1,len(obj))
		if isa(obj.body[i],Matrix)
			CTM=CTM*obj.body[i]
		elseif (isa(obj.body[i],Tuple) || isa(obj.body[i],Array)) && (length(obj.body[i])==2 || length(obj.body[i])==3)
			l=larApply(CTM)(obj.body[i])
			push!(scene,l)
		elseif isa(obj.body[i],Struct)
			push!(stack,CTM)	
			traversal(CTM,stack,obj.body[i],scene)
			CTM=pop!(stack)
			
		end
	end
	return scene
end


#____________________________________________________________________________________________________________________________
function evalStruct(self)
	dim = checkStruct(self.body)
   	CTM, stack = eye(dim+1), []
   	scene = traversal(CTM, stack, self, []) 
return scene
end


