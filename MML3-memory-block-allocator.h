#pragma once
#include<list>
#include<new>

namespace MML3{
/*
La classe MemoryBlockAllocator<T> definisce un allocatore efficiente di memoria specializzato per 
le matrici sparse dinamiche. Si tratta di un allocatore di memoria a blocchi e senza deallocazione atomica
per il tipo template T con le seguenti caratteristiche:
1).	Quando viene richiesta l'allocazione di un elemento, controlla se nel blocco corrente
	c'è spazio. Se si, ritorna il puntatore al primo elemento libero. Se lo spazio è
	esaurito viene allocato un nuovo blocco e viene ritornato il puntatore al primo elemento
	del blocco appena creato. I blocchi esauriti vengono inseriti in una lista.
    La dimensione dei blocchi viene fissata quando l'allocatore viene istanziato e non
	può essere più modificata.
2).	Gli elementi allocati non possono essere singolarmente deallocati. La deallocazione può
	essere eseguita solo sull'intero pool di blocchi.
*/

template<typename T>
class MemoryBlockAllocator
{
public:
	// Costruttore fondamentale a cui si passa la dimensione dei blocchi
	// nei termini del numero di elementi contenuti nel blocco (non in bytes)
	MemoryBlockAllocator( size_t block_size)
		:block_size_(block_size)
	{
		pool_.push_front(new block_(block_size_));
	}
	// Costruttore di default che utilizza blocchi con dimensione 100
	MemoryBlockAllocator()
		:block_size_(100)
	{
		pool_.push_front(new block_(block_size_));
		// std::clog << "allocator contruction 2" << std::endl;
	}

	MemoryBlockAllocator(const MemoryBlockAllocator&) = delete;
	
	// Ritorna il numero di elementi in ciascun blocco
	size_t block_size()const
	{ 
		return block_size_;
	}

	// ritorna un puntatore al primo elemento libero. Se il blocco corrente è esaurito, crea un nuovo blocco.
	// Se il blocco corrente non esiste ancora lo crea.
	T* allocate()
	{
		T* p=pool_.front()->get();
	
		if(!p)
		{
			pool_.push_front(new block_(block_size_));
			p=pool_.front()->get();
		}
		return p;
	}

	~MemoryBlockAllocator()
	{
		std::cout << "allocator destructor called" << std::endl;
		clear();
	}

	// dealloca tutta la memoria allocata
	void clear()
	{
		for (auto& b : pool_)
		{
			b->clear();
			b = nullptr;
		}
		pool_.clear();
		
	}

	size_t	nblocks()const
	{
		return pool_.size();
	}

private:

	struct block_
	{
		block_() = default;
		
		block_(size_t size)
			:begin(new T[size]),size_(size),pos(begin)
		{
		}
		
		~block_(){clear();}
		
		void clear()
		{
			delete [] begin; 
			begin=pos=nullptr;
			size_ = 0;
		}

		T*		get()
		{
			return (pos<begin+size_)?pos++:nullptr;
		}
		size_t	size()const
		{
			returnsize_;
		}

		T*		begin=nullptr;
		size_t  size_=0;
		T*		pos  =nullptr;
	};

	size_t				block_size_=0;	// dimensione del blocco 
	std::list<block_*>	pool_;			// contenitore dei blocchi

}; // end class MemoryBlockAllocator

} // end namespace MML3

