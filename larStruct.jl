#using PyCall
#@pyimport scipy


function prepKey(args)
               v=join(args,",")
return(v)
end


#____________________________________________________________________________________________________________________________

function fixedPrec(PRECISION)
	function fixedPrec0(value)
		out=round(value,PRECISION)
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
		return prepKey(map(fixedPrec(PRECISION),vect))
	end
	return vcode0
end

#____________________________________________________________________________________________________________________________

function t(args...)
	d=length(args)
	mat=eye(d+1)
	for k in range(1,d)
        	mat[k,d+1]=args[k]
	end
	return mat
end

#____________________________________________________________________________________________________________________________

function s(args...)
	d=length(args)
	mat=eye(d+1)
	for k in range(1,d)
		mat[k,k]=args[k]
	end
	return mat
end


#____________________________________________________________________________________________________________________________


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
        mat = eye(3)
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


#____________________________________________________________________________________________________________________________


function larEmbed(k)
	function larEmbed0(model)
		if k>0
			model[1]=[append!(v,fill(0.0,k)) for v in model[1]]
		elseif k<0
			for n in range(1,length(model[1][1])+k)
				for v in model[1]
					pop!(v)
				end
			end
		end
		return model
	end
	return larEmbed0
end


#____________________________________________________________________________________________________________________________

function removeDups(CW)
	CW=collect(tuple(CW...))
	CWs=collect(tuple(sort(CW,rev=true)...))
	no_duplicates=Dict()
	for f in CWs
		no_duplicates[f] = []
	end
	for f in CW
	no_duplicates[tuple(sort(collect(f))...)]=[f]
	end
	CW=[f[1] for f in values(no_duplicates)]
	return CW
end 


#____________________________________________________________________________________________________________________________

 
function larRemoveVertices(V,FV)
	vertDict= Dict()
	index,defaultValue,CW,W,FW = 0,-1,[],[],[]
	for (k,incell) in enumerate(FV)
		outcell=[]
		for v in incell
			key=vcode(4)(V[v])
			if get(vertDict,key,defaultValue)==defaultValue
				index =index+1
				vertDict[key]=index
				outcell =append!(outcell,index)
				W= append!(W,eval(key))                   
			else
				outcell= append!(outcell,vertDict[key])
			end
		end
		FW =append!(FW,[outcell])
	end
	return W,FW
end



#____________________________________________________________________________________________________________________________


#classe Struct



type Struct
	body::Array
	box::Array
	name::AbstractString
	dim::Int
	category::AbstractString
	Struct()=new([],Nullable,string(object_id(new())),Nullable,"feature")
	Struct(data::Array)=new([item for item in data],box(new(),data),string(object_id(new())),length(box[1]),"feature")
	Struct(data::Array,name)=new([item for item in data],box(this),string(name),length(box[1]),"feature")
	Struct(data::Array,name,category)=new([item for item in data],box(this),string(name),length(box[1]),string(category))
end

	function __name__(self::Struct)
		return self.name
	end
	function __category__(self::Struct)
		return self.category
	end
	#function __iter__(self::Struct)
		#return take(self.body,length(self.body))
	#end
	function __len__(self::Struct)
		return length(self.body)
	end
	function __getitem__(self::Struct,i::Int)
		return self.body[i]
	end
	function __setitem__(self::Struct,i,value)
		self.body[i]=value
	end
	function __print__(self::Struct)
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

function larBoundary(self)
	data=struct2lar(self)
	if length(data)==3
		V,FV,EV=data
		return V,FV,EV
	else
		return "<Struct name $(self.__name__()) : boundary non computable"
	end
end


#____________________________________________________________________________________________________________________________

function struct2lar(structure,metric=ID)
	listOfModels=evalStruct(structure)
	vertDict= Dict()
	index,defaultValue,CW,W,FW = 0,-1,[],[],[]
	for model in listOfModels
		if  length(model)==2
			V,FV=model
		elseif lenght(model)==3
			V,FV,EV=model
		end
		for (k,incell) in enumerate(FV)
			outcell=[]
			for v in incell
				key=vcode(4)(V[v])
				if get(vertDict,key,defaultValue)==defaultValue
					index =index+1
                                        vertDict[key]=index
					outcell =append!(outcell,index)
					W= append!(W,eval(key))                   
				else
					outcell= append!(outcell,vertDict[key])
				end
			end
			CW =append!(CW,[outcell])
		end
		if length(model)==3
			for (k,incell) in enumerate(FV)
				outcell=[]
				for v in incell
					key=vcode(4)(V[v])
					if get(vertDict,key,defaultValue)==defaultValue
						index =index+1
						vertDict[key]=index
						outcell =append!(outcell,[index])
						W=append!(W,eval(key))                   
					else
						outcell= append!(outcell,vertDict[key])
					end
				end
				FW =append!(FW,[outcell])
			end
		end
	end
	if length(model)==2
		elseif length(CW[1])==2
			CW=map(Tuple,map(sort,CW))
		else
			CW=removeDups(CW)
		return metric(W),CW
	end
	if length(model)==3
		FW=map(Tuple,map(sort,FW))
		CW=removeDups(CW)
		return metric(W),CW,FW
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
			newObj.box=hcat((obj[i].box,[fill([0],n),fill([0],n)]))	#vedere con test hosppital2/01 se viene lo stesso risultato
			newObj.category=obj[i].category
			append!(cloned.body,embedTraversal(newObj,obj[i],n,suffix))
		end
	end
	return cloned
end

#____________________________________________________________________________________________________________________________


function embedStruct(n)
	function embedStruct0(struct,suffix="New")
		if n==0
			return struct, length(struct.box[1])
		end
		cloned=Struct()
		cloned.box=hcat((struct.box,[fill([0],n),fill([0],n)]))	#vedere con test hosppital2/01 se viene lo stesso risultatok
		cloned.name=string(object_id(cloned))
		cloned.category=struct.category
		cloned.dim=struct.dim+n
		cloned=embedTraversal(cloned,struct,n,suffix)
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
		dummyModel.body=[]
		for term in model.body
			if isa(term,Struct)
				append!(dummyModel.body,[term.box,Array[0,1]])	#se da errore provare: ,Array[term.box,ecc...]
			else
				append!(dummyModel.body,term)
			end
		end
		listOfModels=evalStruct(dummyModel)
		#dim=checkStruct(listOfModels)
		theMin,theMax=box(listOfModels[0])
		for theModel in listOfModels[2:end]
			modelMin,modelMax= box(theModel)
			for (k,val) in enumerate(modelMin)
				if val<theMin[k]
					theMin=[val]
				else
					theMin=theMin[k]
				end
			end
			for (k,val) in enumerate(modelMax)
				if val>theMax[k]
					theMax=[val]
				else
					theMax=theMax[k]
				end
			end
		end
		return[theMin,theMax]
	elseif (isa(model,Tuple) ||isa(model,Array))&& (length(model)==2 || length(model)==3)
		V=model[1]
	end
	theMin=[minimum(V[1:end,i]) for i in range(1,ndims(V)+1)]
	theMax=[maximum(V[1:end,i]) for i in range(1,ndims(V)+1)]
	return Array[theMin,theMax]
end

 

 

#____________________________________________________________________________________________________________________________


function larApply(affineMatrix)
  function larApply0(model)
    if length(model)==2
      V,CV=model
    elseif length(model)==3
      V,CV,FV = model
    end
	 for (k,v) in enumerate(V)
		append!(v,[1])
    		V[k]=vec((v')*(transpose(affineMatrix)))
	end

  if length(model)==2
	for v in V
		pop!(v)
	end
	return V,CV	#una cosa simile la avevo anche io...ho fatto un for(in larEmbed),ricordiamoci di confrontare!

  elseif length(model)==3
	V,CV,FV = model
	for v in V
		pop!(v)
	end
	return V,CV,FV

  end

end 

return larApply0

end

#____________________________________________________________________________________________________________________________

function checkStruct(lst)
	obj = lst[1]
	if(isa(obj,Tuple) || isa(obj,Array))
		dim=length(obj[1][1])
	elseif isa(obj,Matrix)
		dim=size(obj)[1]-1
	elseif isa(obj,Struct)
		dim=length(obj.box[1])
	end
	return dim
end		





#____________________________________________________________________________________________________________________________


function traversal(CTM,stack,obj,scene=[])
	for i in range(1,length(obj))
		if((isa(obj[i],Tuple) || isa(obj[i],Array)) & (length(obj[i]==2) || length(obj[i]==3)))
			apend!(scene,larApply(CTM)(obj[i]))
		elseif(isa(obj[i],Matrix))
			CTM=CTM*obj[i]
		elseif(isa(obj[i],Struct))
			append!(CTM,stack)	#non va bene dobbiamo trovare il modo per mettere in ogni posizione di stack una matrice. 
			traversal(CTM,stack,obj,scene)
			#stessa cosa per togliere l'ultimo elemento di stack, che però è una matrice CTM=stack.pop()
		end
	end
	return scene
end


#____________________________________________________________________________________________________________________________
function evalStruct(struct)
	dim = checkStruct(struct.body)
   	CTM, stack = eye(dim+1), []
   	scene = traversal(CTM, stack, struct, []) 
return scene
end












	
