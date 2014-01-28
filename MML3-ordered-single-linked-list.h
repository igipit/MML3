#pragma once

#include<iomanip>

namespace MML3
{


// forward declarations
template<typename VAL, typename IDX>
class list_iterator;

template<typename VAL, typename IDX>
class const_list_iterator;

/*
 La componente di una lista concatenata semplice
 Con allineamento a 8 bytes, architettura 32 bit, sparse_component<double,__int32> occupa 
 16 bytes che è il minimo possibile
*/


template<typename VAL, typename IDX>
struct sparse_component
{
	typedef VAL										value_t;
	typedef IDX										index_t;
	typedef list_iterator<value_t,index_t>			iterator;
	typedef const_list_iterator<value_t,index_t>	const_iterator;

			sparse_component(){	}
			
			sparse_component(value_t	val, index_t i, sparse_component*	nxt = nullptr) :value(val), idx(i), next(nxt){}

			sparse_component* set(value_t	val, index_t i, sparse_component*	nxt = nullptr)
			{
				value=val;
				idx=i;
				next=nxt;
				return this;
			}
			
	value_t					value=0;
	index_t					idx=0;	
	sparse_component		*next = nullptr;
};

template<typename VAL, typename IDX>
class list_iterator
{
	friend class const_list_iterator<VAL,IDX>;

	public:
		typedef sparse_component<VAL,IDX> item;
		list_iterator(item* init):p_(init){}
		list_iterator(){}

		
		list_iterator&	operator ++(){p_=p_->next;return *this;}
		list_iterator	operator ++(int){list_iterator pro(p_); p_=p_->next; return pro;}
		list_iterator&	operator =(const list_iterator& rhs){p_=rhs.p_; return *this;}
		item*			operator->(){return p_;}
		item&			operator*(){return *p_;}
		bool			operator!=(const list_iterator& rhs)const{return p_!=rhs.p_;}
		bool			operator!=(const const_list_iterator<VAL,IDX>& rhs)const{return p_!=rhs.p_;}
		bool			is_valid()const{ return p_ != nullptr; }
		operator bool(){ return p_ != nullptr; }
		VAL&			value(){ return p_->value; }
		const IDX&		idx()const{ return p_->idx;}

	private:
		item* p_ = nullptr;
};

template<typename VAL, typename IDX>
inline bool	is_valid(const list_iterator<VAL,IDX>& lit){return lit.is_valid();}



template<typename VAL, typename IDX>
class const_list_iterator
{
	friend class list_iterator<VAL,IDX>;
	public:
		typedef sparse_component<VAL,IDX> item;
		const_list_iterator(const item* init):p_(init){}
		const_list_iterator(){}
		const_list_iterator(const list_iterator<VAL,IDX>& rhs):p_(rhs.p_){}
		
		
		const_list_iterator&	operator ++(){p_=p_->next;return *this;}
		const_list_iterator	operator ++(int){list_iterator<VAL,IDX> pro(p_); p_=p_->next; return pro;}
		const_list_iterator&	operator =(const list_iterator<VAL,IDX>& rhs){p_=rhs.p_; return *this;}
		const_list_iterator&	operator =(const const_list_iterator& rhs){p_=rhs.p_; return *this;}
		const item*				operator->()const{return p_;}
		const item&				operator*()const{return *p_;}
		bool					operator!=(const const_list_iterator& rhs)const{return p_!=rhs.p_;}
		bool					is_valid()const{ return p_ != nullptr; }
		const VAL&				value(){ return p_->value; }
		const IDX&				idx()const{ return p_->idx; }
		operator bool(){ return p_ != nullptr; }
	private:
		mutable item* p_ = nullptr;
};

template<typename VAL, typename IDX>
inline bool	is_valid(const const_list_iterator<VAL,IDX>& lit){return lit.is_valid();}




/* 
Lista concatenata semplice ordinata di elementi  sparse_component<value_t,index_t>
*/

template<typename VAL, typename IDX, class ALLOCATOR >
class ordered_single_linked_list
{
public:

	typedef sparse_component<VAL,IDX>	element;
	typedef element*					element_pointer;
	typedef element_pointer*			element_pointer_pointer;
	typedef VAL							value_t;
	typedef IDX							index_t;
	typedef typename element::const_iterator	const_iterator;
	typedef typename element::iterator	iterator;
	

	ordered_single_linked_list(ALLOCATOR* alloc):head_(),allocator_(alloc){}
	ordered_single_linked_list(ALLOCATOR* alloc, const ordered_single_linked_list& rhs);

	// !!! Attenzione, la lista non dealloca mai la memoria allocata
	~ordered_single_linked_list(){ }

	// conta e ritorna il numero di elementi nella lista operatore O(1)
	index_t	size()const{return head_.idx;}
	

	/////////////////////////////////
	// COMPONENT SEARCH & ACCESS
	iterator		begin(){return iterator(head_.next);}
	const_iterator	begin()const{return iterator(head_.next);}
	const_iterator	end()const{return const_iterator();}
	iterator		end(){ return iterator(); }

	element*		first(){return head_.next;}
	const element*	first()const{return head_.next;}

	
	// Cerca l'elemento precedente a quello di indice idx. Se quest'ultimo non c'è torna il puntatore all'ultimo elemento della lista.
	// Questa funzione ritorna sempre un puntatore valido ad un elemento della lista
	element* find_first_before(index_t idx)
	{
		element* p=&head_;
		for(; p->next && p->next->idx<idx; p=p->next);
		return p;
	}

	// versione const della precedente
	const element* find_first_before(index_t idx)const
	{
		const element* p=&head_;
		for(; p->next && p->next->idx<idx; p=p->next);
		return p;
	}


	// ritorna il puntatore all'elemento di indice idx. Se questo non c'è torna nullptr
	const element*	get(index_t idx)const
	{
		const element* p=find_first_before(idx);
		if(p->next->idx==idx)
			return p->next;
		else
			return nullptr;
	}
	
	element*	get(index_t idx)
	{
		const element* p = find_first_before(idx);
		if (p->next->idx == idx)
			return p->next;
		else
			return nullptr;
	}
	/*
	//cerca l'elemento di indice i. Se lo trova inserisce il valore in val e pone found=true. Altrimenti found=false
	void	get(index_t i, value_t& val, bool& found)const
	{
		const element*	p=get(i);
		if(p)
		{
			val=p->value;
			found=true;
		}
		else
			found =false;

	}
	*/

	// gets the first n element of the list and put them into the i_buffer and v_buffer
	// return the number of element inserted that is min(n, size())
	int	get_first_n(index_t* i_buffer , value_t* v_buffer, index_t n)const
	{
		size_t m=std::min(n,size());
		const element* p=head_.next;;
		for(size_t i=0; i!=m; ++i)
		{
			if(i_buffer)
				i_buffer[i]=p->idx;
			if (v_buffer)
				v_buffer[i]=p->value;
			p=p->next;
		}
		return m;
	}

	void	get(std::vector<index_t>& i_vec,std::vector<value_t>& v_vec)const
	{
		size_t n=head_.idx;
		i_vec.resize(n);
		v_vec.resize(n);
		if (n)
			get_first_n(&i_vec[0], &v_vec[0], n);
		
	}


	void	get_indexes(std::vector<index_t>& i_vec)const
	{
		size_t n=head_.idx;
		i_vec.resize(n);
		if (n)
			get_first_n(&i_vec[0], nullptr, n);
	}

	void	get_values(std::vector<value_t>& v_vec)const
	{
		size_t n=head_.idx;
		v_vec.resize(n);
		if (n)
			get_first_n(nullptr, &v_vec[0], n);
		
	}

	
	void	fill(value_t val)
	{
		for( element* p = head_.next; p; p=p->next)
			p->value=val;
	}
	
	

	// Inserisce nella lista un  elemento di indice idx. Se questo non c'è lo crea
	// altrimenti sovrascrive il valore gia' presente
	value_t* insert(index_t idx, const value_t& val)
	{
		element* p=find_first_before(idx);
		if(!p->next || p->next->idx!=idx)
		{
			p->next=allocator_->allocate()->set(val,idx,p->next);
			++head_.idx;
		}
		else
			p->next->value=val;
		
		return &(p->next->value);

	}

	

	// Aggiunge all'elemento di indice idx il valore val. Se l'elemento non c'è lo crea
	value_t* add_at(index_t idx, const value_t& val)
	{
		element* p=find_first_before(idx);
		if(!p->next || p->next->idx!=idx)
		{
			p->next=allocator_->allocate()->set(val,idx,p->next);
			++head_.idx;
		}
		else
			p->next->value+=val;

		return &(p->next->value);
	}

	
	// riempie la lista con n elementi (nell'ordine  in cui sono dati)
	// attenzione, il contenuto precedente viene perso e le relative risorse non deallocate
	void set(const index_t* i_buff, const value_t* v_buff, size_t sz)
	{
		element* p=&head_;
		head_.idx = 0;
		for(size_t i=0; i!=sz; ++i)
		{
			if (!p->next)
				p->next = allocator_->allocate()->set(v_buff[i], i_buff[i], nullptr);
			else
				p->next->set(v_buff[i], i_buff[i], p->next->next);
			++head_.idx;
			p=p->next;
		}
		p->next = nullptr;
	}
	


	///////////////////////////////////////////////////////////////////////////////
	// IO
	// Stampa in output la lista. Sulla prima riga gli indici, sulla seconda i valori
	
	void print(std::ostream& os, unsigned int width)
	{
	
		int prec=os.precision();
		for(element* p=head_.next; p; p=p->next)
			os	<< std::setw(width) << (unsigned int) p->idx ;
		os	<< "\n";
		for(element* p=head_.next; p; p=p->next)
			os	<< std::setw(width) <<  std::setprecision(prec) << p->value ;
	}

	

	
	// inserisce nella testa della lista il contenuto di un buffer di valori ed uno di indici
	// questi devono essere già ordinati altrimenti la lista risulterà disordinata. Inoltre la lista deve essere inzialmente vuota
	// altrimenti viene lanciata un'eccezione
	/*
	void insert(value_t* v_buff, index_t* i_buff, index_t size)
	{
		element*	p=head_.nxt;
		if(v_buff && i_buff)
			for(index_t i=size; i; --i)
				p=allocator_->allocate()->set(v_buff[i-1],i_buff[i-1],p);
		else if(v_buff)
				for(index_t i=size; i; --i)
					p=allocator_->allocate()->set(v_buff[i-1],0,p);
		else if(i_buff)
				for(index_t i=size; i; --i)
					p=allocator_->allocate()->set(0,i_buff[i-1],p);
	}
	*/
private:

	// head è un elemento artificiale il cui puntartore next punta al primo elemento della lista. Il campo idx viene usato per contare gli elementi della lista
	element		head_;
	ALLOCATOR*	allocator_=nullptr;
};


template<typename VAL, typename IDX, class ALLOCATOR >
ordered_single_linked_list<VAL,IDX,ALLOCATOR>::ordered_single_linked_list(ALLOCATOR* alloc, const ordered_single_linked_list& rhs) :head_(), allocator_(alloc)
{
	element			*dest = &head_;
	for (const element *src = rhs.head_.next; src; src = src->next)
	{
		dest->next = allocator_->allocate()->set(src->value, src->idx, nullptr);
		dest = dest->next;
	}
}

} // end namespace MML3
