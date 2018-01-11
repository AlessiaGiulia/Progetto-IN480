using PyCall
@pyimport scipy


function prepKey(args)
               v=join(args,",")
return(v)
end


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

 

function vcode(PRECISION=4)
	function vcode0(vect)
		return prepKey(map(fixedPrec(PRECISION),vect))
	end
	return vcode0
end


function t(args…)
	d=length(args)
	mat=eye(d+1)
	for k in range(1,d)
        	mat[k,d+1]=args[k]
	end
	return mat
end

 

function s(args…)
	d=length(args)
	mat=eye(d+1)
	for k in range(1,d)
		mat[k,k]=args[k]
	end
	return mat
end


function r(args)
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
        if axis[2]==axis[3]==0.0:    # rotation about x
            mat[2,2] = COS;    mat[2,3] = -SIN;
            mat[3,2] = SIN;    mat[3,3] = COS;

        elseif axis[1]==axis[3]==0.0:    # rotation about y
            mat[1,1] = COS;    mat[1,3] = SIN;
            mat[3,1] = -SIN;    mat[3,3] = COS;
        elseif axis[0]==axis[1]==0.0:    # rotation about z
            mat[0,0] = SIN;    mat[0,1] = -SIN;
            mat[1,0] = COS;    mat[1,1] = COS;
        
        else
	    I=eye(3); u=axis
	    Ux=[[0,-u[3],u[2]],  [u[3],0,-u[1]],  [-u[2],u[1],1]]
	    UU =[[u[1]*u[1],    u[1]*u[2],    u[1]*u[3]],
                 [u[2]*u[1],    u[2]*u[2],    u[2]*u[3]],
                 [u[3]*u[1],    u[3]*u[2],    u[3]*u[3]]])
	    mat[1:3,1:3]=COS*I+SIN*Ux+(1.0-COS)*UU
	end

return mat



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
				outcell =append!(outcell,index])
				W= append!(W,eval(key))                   
			else
				outcell= append!(outcell,vertDict[key])
			end
		end
		FW =append!(FW,[outcell])
	end
	return W,FW
end




#classe Struct
type Struct
	body::Array
	box::Array
	name::String
	dim:: Int
	category::String
	
	function Struct(self,data=None,name=None,category=None)
		self=new()
		if (data==None || data=[])
			self.body=[]
		else
			self.body=[item for item in data]
			self.box=box(self)
			self.dim=length(self.box[1])
		end
		if name!= None
			self.name=string(name)
		else
			self.name=string(object_id(self))
		end
		if category!= None
			self.category=string(category)
		else
			self.category= "feature"
		end
	end

	function __name__(self)
		return self.name
	end
	function __category__(self)
		return self.category
	end
	#function __iter__(self)
		#return take(self.body,length(self.body))
	#end
	function __len__(self)
		return length(self.body)
	end
	function __getitem__(self,i)
		return self.body[i]
	end
	function __setitem__(self,i,value)
		self.body[i]=value
	end
	function __print__(self)
		return "<Struct name: $(self.__name__())"
	end
	function set_name(self,name)
		self.name=string(name)
	end
	function clone(self,i=0)
		newObj=deepcopy(self)
		if i!=0
			newObj.name="$(self.__name__())_$(string(i))"
		end
		return newObj
	end
	function set_category(self,category)
		self.category=string(category)
	end
end


function larBoundary(self)
	data=struct2lar(self)
	if length(data)==3
		V,FV,EV=data
		return V,FV,EV
	else
		return "<Struct name $(self.__name__()) : boundary non computable"
	end
end


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
					outcell =append!(outcell,index])
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
		elseif length(CW[0])==2
			CW=list(set(AA(tuple)(AA(sorted)(CW))))#da cambiare
		else
			CW=removeDups(CW)
		return metric(W),CW
	end
	if length(model)==3
		FW=list(set(AA(tuple)(AA(sorted)(FW))))#da cambiare
		CW=removeDups(CW)
		return metric(W),CW,FW
	end
end

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


function box(model)
	if isa(model,Matrix)
		return []
	elseif isa(model,Struct)
		dummyModel=deepcopy(model)
		dummyModel.body=[]
		for term in model.body
			if isa(term,Struct)
				append!(dummyModel.body,[term.box,Array[0,1]]	#se da errore provare: ,Array[term.box,ecc...]
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
	elseif (isa(model,tuple) ||isa(model,Array))&& (length(model)==2 || length(model)==3)
		V=model[1]
	end
	coord=hcat(V...)
	coord=[coord[i,:] for i in range(1,length(V[1]))]
	theMin=min(coord...)
	theMax=max(coord...)
	return Array[theMin,theMax]
end

 

 

--------------------------------------------------------------------------------------------------------------------------------------------

function larApply(affineMatrix)
  function larApply0(model)
    if length(model)==2
      V,CV=model
    elseif length(model)==3
      V,CV,FV = model
    end

    for v in V 
      append!(v,1.0) 
    end  
 
    V=(v')*(transpose(affineMatrix))

  if len(model)==2
	for v in V
		pop!(v)
	end
	return V,CV	#una cosa simile la avevo anche io...ho fatto un for(in larEmbed),ricordiamoci di confrontare!

  elseif len(model)==3
	for v in V
		pop!(v)
	end
	return V,CV,FV

  end

end 

return larApply0

end

#___________________________________________________________________________

function checkStruct(lst)
	obj = lst[1]
	





#_________________________________________________________________________


function traversal(CTM,stack,obj,scene=[])




#______________________________________________________________________
function evalStruct(struct)
	dim = checkStruct(struct.body)
   	CTM, stack = eye(dim+1), []
   	scene = traversal(CTM, stack, struct, []) 
return scene
end













	
