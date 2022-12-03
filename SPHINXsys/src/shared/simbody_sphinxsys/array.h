/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4													*
 *                                                                          *
 * Portions copyright (c) 2017-2022 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	array.h
 * @details An Array toolkit for array operator. 
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef ARRAY_H
#define ARRAY_H

#include "base_data_package.h"

#include <iostream>
#include <string>
#include <cstdio>
#include <stdexcept>

static const int Array_CAPMIN = 1;

namespace SPH {
	/**
	 * @class Array
 	 * @ details A class for storing an array of values of type T.  The capacity of the class
 	 * 		grows as needed.  To use this template for a class of type T, class T should
 	 * 		implement the following methods:  default constructor, copy constructor,
 	 * 		assignment operator (=), equality operator (==), and less than
 	 * 		operator (<).
 	 */
	template<class T> class Array
	{
	protected:
    	/** Size of the array.  Also the index of the first empty array element. */
    	int size_;
    	/** Current capacity of the array. */
    	int capacity_;
    	/** Increment by which the current capacity is increased when the capacity
    	of this storage instance is reached.  If negative, capacity doubles. */
    	int capacityIncrement_;
    	/** Default value of elements. */
    	T defaultValue_;
    	/** Array of values. */
    	T *array_;
    public:
		/**
 		* @brief Default constructor.
 		*
 		* @param aDefaultValue Default value of an array element.  This value
 		* 		is used to initialize array elements as the size of the array is
 		* 		changed.
 		* @param aSize Initial size of the array.  The array elements are
 		* 		initialized to aDefaultValue.
 		* @param aCapacity Initial capacity of the array.  The initial capacity
 		* 		is guaranteed to be at least as large as aSize + 1.
 		*/
		explicit Array(const T &aDefaultValue=T(),int aSize=0,int aCapacity=Array_CAPMIN)
		{
    		setNull();
    		defaultValue_ = aDefaultValue;
    		int newCapacity;
    		int min = aSize + 1;
    		if(min < aCapacity) min = aCapacity;
    		computeNewCapacity(min,newCapacity);
    		ensureCapacity(newCapacity);
    		size_ = aSize;
    		if(size_<0) size_=0;
		}
		/**
 		 * @brief Copy constructor.
 		 *
 		 * @param aArray Array to be copied.
 		 */
		Array(const Array<T> &aArray)
		{
    		setNull();
    		*this = aArray;
		}
		/**
 		 * @brief 	Destructor.
 		 *
 		 * @details When the array is deleted, references to elements of this array become
 		 * 		  	invalid.
 		 */
		virtual ~Array()
		{
    		if(array_!=nullptr) { delete[] array_;  array_ = nullptr; }
		}
	private:
		/**
 		* %Set all member variables to their null or default values.
 		*/
		void setNull()
		{
    		size_ = 0;
    		capacityIncrement_ = -1;
    		capacity_ = 0;
    		array_ = nullptr;
		}
	//=============================================================================
	// OPERATORS
	//=============================================================================
	public:
		/**
		 * @brief A non-operator version of operator == 
		 */
		bool arrayEquals(const Array<T> &aArray) const
		{
    		return *this == aArray;
		}
		/**
 		 * @brief Get the array element at a specified index.  This overloaded operator
 		 * 			can be used both to set and get element values:
 		 * @code
 		 *      Array<T> array(2);
 		 *      T value = array[i];
 		 *      array[i] = value;
 		 * @endcode
 		 *
 		 * @details This operator is intended for accessing array elements with as little
 		 * 			overhead as possible, so no error checking is performed.
 		 * 			The caller must make sure the specified index is within the bounds of
 		 * 			the array.  If error checking is desired, use Array::get().
 		 *
 		 * @param aIndex Index of the desired element (0 <= aIndex < size_).
 		 * @return Reference to the array element.
 		 * @see get().
 		*/
		T& operator[](int aIndex) const
		{
    		return(array_[aIndex]);
		}
		/**
 		 * @brief Assign this array to a specified array.
 		 * 		This operator makes a complete copy of the specified array; all member
 		 * 		variables are copied.  So, the result is two identical, independent arrays.
 		 *
 		 * @param aArray Array to be copied.
 		 * @return Reference to this array.
 		*/
		Array<T>& operator=(const Array<T> &aArray)
		{
    		size_ = aArray.size_;
    		capacity_ = aArray.capacity_;
    		capacityIncrement_ = aArray.capacityIncrement_;
    		defaultValue_ = aArray.defaultValue_;

    		// ARRAY
    		if(array_!=nullptr) delete[] array_;
    		array_ = new T[capacity_];
    		for(int i = 0;i < capacity_; i++) array_[i] = aArray.array_[i];

    		return(*this);
		}
		/**
 		 * @brief Determine if two arrays are equal.
 		 *
 		 * @details Two arrays are equal if their contents are equal.  That is, each array
 		 * 			must be the same length and their corresponding array elements must be
 		 * 			equal.
 		 *
 		 * @param aArray Array to be tested as equal.
 		 * @return True if equal, false if not equal.
 		*/
		bool operator==(const Array<T> &aArray) const
		{
    		if(size_ != aArray.size_) return(false);

    		int i;
    		for(i=0; i<size_; i++) 
    		{
        		if( !(array_[i]==aArray.array_[i]) ) return(false);
    		}

    		return(true);
		}
		/**
 		 * @brief Implementation of the output operator.
 		 * 			The output for an array looks like the following:\n\n
 		 *
 		 * T[0] T[1] T[2] ... T[size-1].
 		 *
 		 * @param aOut Output stream.
 		 * @param aArray Array to be output.
 		 * @return Reference to the output stream.
 		*/
		friend std::ostream& operator<<(std::ostream &aOut,const Array<T> &aArray)
		{
    		int i;
    		for(i=0; i<aArray.getSize(); i++)  {
        		aOut << " ";
        		aOut << aArray[i];
    		}	

    		return(aOut);
		}

		friend std::istream& operator>>(std::istream& in, Array<T>& out) 
		{    
    		return in;
		}
		/**
 		 * @brief Compute a new capacity that is at least as large as a specified minimum
 		 * 			capacity; this method does not change the capacity, it simply computes
 		 * 			a new recommended capacity.
 		 *
 		 * @details If the capacity increment is negative, the current capacity is
 		 * 			doubled until the computed capacity is greater than or equal to the
 		 * 			specified minimum capacity.  If the capacity increment is positive, the
 		 * 			current capacity increment by this amount until the computed capacity is
 		 * 			greater than or equal to the specified minimum capacity.  If the capacity
 		 * 			increment is zero, the computed capacity is set to the current capacity
 		 * 			and false is returned.
 		 *
 		 * @param rNewCapacity New computed capacity.
 		 * @param aMinCapacity Minimum new computed capacity.  The computed capacity
 		 * 			is incremented until it is at least as large as aMinCapacity, assuming
 		 * 			the capacity increment is not zero.
 		 * @return True if the new capacity was increased, false otherwise (i.e.,
 		 * 			if the capacity increment is set to 0).
 		 * @see setCapacityIncrement()
 		*/
		bool computeNewCapacity(int aMinCapacity,int &rNewCapacity)
		{
    		rNewCapacity = capacity_;
    		if(rNewCapacity < Array_CAPMIN) rNewCapacity = Array_CAPMIN;

    		if(capacityIncrement_ == 0) {
        		std::cout << "Array.computeNewCapacity: WARN- capacity is set";
        		std::cout << " not to increase (i.e., capacityIncrement_==0).\n";
        		return(false);
    		}

    		while(rNewCapacity < aMinCapacity) {
        		if(capacityIncrement_ < 0) {
            		rNewCapacity = 2 * rNewCapacity;
        		} else {
            		rNewCapacity = rNewCapacity + capacityIncrement_;
        		}
    		}

    		return(true);
		}
		/**
 		* @brief Ensure that the capacity of this array is at least the specified amount.
 		* 		Note that the newly allocated array elements are not initialized.
 		*
 		* @param aCapacity Desired capacity.
 		* @return true if the capacity was successfully obtained, false otherwise.
 		*/
		bool ensureCapacity(int aCapacity)
		{
    		if(aCapacity < Array_CAPMIN) aCapacity = Array_CAPMIN;
    		if(capacity_ >= aCapacity) return(true);

    		int i;
    		T *newArray = new T[aCapacity];
    		if(newArray == nullptr) 
    		{
        		std::cout << "Array.ensureCapacity: ERR- failed to increase capacity.\n";
        		return(false);
    		}

    		if(array_ != nullptr) {
        		for(i =0; i < size_; i++) newArray[i] = array_[i];
        		for(i = size_; i < aCapacity; i++) newArray[i] = defaultValue_;
        		delete []array_;  
        		array_ = nullptr;
    		} else {
        		for(i = 0; i < aCapacity; i++) newArray[i] = defaultValue_;
    		}
    
    		capacity_ = aCapacity;
    		array_ = newArray;

    		return(true);
		}

		/**
 		* @details Trim the capacity of this array so that it is one larger than the size
 		* 	of this array.  This is useful for reducing the amount of memory used
 		* 	by this array.  This capacity is kept at one larger than the size so
 		* 	that, for example, an array of characters can be treated as a nullptr
 		* 	terminated string.
 		*/
		void trim()
		{
    		int newCapacity = size_ + 1;
    		if(newCapacity >= capacity_) return;
    		if(newCapacity < Array_CAPMIN) newCapacity = Array_CAPMIN;

    		int i;
    		T *newArray = new T[newCapacity];
    		if(newArray==nullptr) {
        		std::cout << "Array.trim: ERR- unable to allocate temporary array.\n";
        		return;
    		}

    		for(i = 0; i < size_; i++) newArray[i] = array_[i];

    		delete[] array_;

    		array_ = newArray;

    		capacity_ = newCapacity;
		}
		/**
 		 * Get the capacity of this storage instance.
 		 */
		int getCapacity() const
		{

    		return(capacity_);
		}
		/**
 		* @brief Set the amount by which the capacity is increased when the capacity of
 		* 		of the array in exceeded.
 		* @details If the specified increment is negative, the capacity is set to double
 		* 			whenever the capacity is exceeded.
 		*
 		* @param aIncrement Desired capacity increment.
 		*/
		void setCapacityIncrement(int aIncrement)
		{
    		capacityIncrement_ = aIncrement;
		}

		/**
 		 * @brief Get the amount by which the capacity is increased.
 		*/
		int getCapacityIncrement() const
		{
    		return(capacityIncrement_);
		}
		/**
 		 * @details Set the size of the array.  This method can be used to either increase
 		 * 			or decrease the size of the array.  If this size of the array is
 		 * 			increased, the new elements are initialized to the default value
 		 * 			that was specified at the time of construction.
 		 *
 		 * @details Note that the size of an array is different than its capacity.  The size
 		 * 			indicates how many valid elements are stored in an array.  The capacity
 		 * 			indicates how much the size of the array can be increased without
 		 * 			allocated more memory.  At all times size <= capacity.
 		 *
 		 * @param aSize Desired size of the array.  The size must be greater than
 		 * or equal to zero.
 		*/
		bool setSize(int aSize)
		{
    		if(aSize == size_) return(true);
    		if(aSize < 0) aSize = 0;
    		bool success = true;
    		if(aSize < size_) 
    		{
        		int i;
        		for(i = (size_ - 1);i >= aSize; i--) array_[i] = defaultValue_;
        		size_ = aSize;
    		} else if(aSize <= capacity_) {
        		size_ = aSize;
    		} else {
        		int newCapacity;
        		success = computeNewCapacity(aSize+1, newCapacity);
        		if(!success) return(false);
        		success = ensureCapacity(newCapacity);
        		if(success) size_ = aSize;
    		}

    		return(success);
		}
		/**
 		 * @brief Get the size of the array.
 		 *
 		 * @return Size of the array.
 		*/
		int getSize() const
		{
    		return(size_);
		}
		/** Alternate name for getSize(). **/
		int size() const {return getSize();}

		/**
 		 * @brief Append a value onto the array.
 		 *
 		 * @param aValue Value to be appended.
 		 * @return New size of the array, or, equivalently, the index to the new
 		 * 			first empty element of the array.
 		*/
		int append(const T &aValue)
		{
    		if((size_ + 1) >= capacity_) {
        		int newCapacity;
        		bool success;
        		success = computeNewCapacity(size_ + 1, newCapacity);
        		if(!success) return(size_);
        		success = ensureCapacity(newCapacity);
        		if(!success) return(size_);
    		}

    		array_[size_] = aValue;
    		size_++;

    		return(size_);
		}
		/**
 		 * @brief Append an array of values.
 		 *
 		 * @param aArray Array of values to append.
 		 * @return New size of the array, or, equivalently, the index to the new
 		 * 			first empty element of the array.
 		*/
		int append(const Array<T> &aArray)
		{
    		int i,n = aArray.getSize();
    		for(i = 0; i < n; i++) {
        		append(aArray[i]);
    		}

    		return(size_);
		}

		/**
 		 * @brief 	Append an array of values.
 		 *
 		 * @param aSize Size of the array to append.
 		 * @param aArray Array of values to append.
 		 * @return New size of the array, or, equivalently, the index to the new
 		 * first empty element of the array.
 		*/
		int append(int aSize,const T *aArray)
		{
    		if(aSize < 0) return(size_);
    		if(aArray == nullptr) return(size_);

    		int i;
    		for(i = 0;i < aSize; i++) {
        		append(aArray[i]);
    		}

    		return(size_);
		}
		/**
 		 * @brief Insert a value into the array at a specified index.
 		 *
 		 * @details This method is relatively computationally costly since many of the array
 		 * 			elements may need to be shifted.
 		 *
 		 * @param aValue Value to be inserted.
 		 * @param aIndex Index at which to insert the new value.  All current elements
 		 * 		from aIndex to the end of the array are shifted one place in the direction
 		 * 		of the end of the array.  If the specified index is greater than the
 		 * 		current size of the array, the size of the array is increased to aIndex+1
 		 * 		and the intervening new elements are initialized to the default value that
 		 * 		was specified at the time of construction.
 		 * @return Size of the array after the insertion.
 		*/
		int insert(int aIndex,const T &aValue)
		{
    		if(aIndex<0) {
        		std::cout << "Array.insert: ERR- aIndex was less than 0.\n";
        		return(size_);
    		}

    		if(aIndex >= size_) {
        		setSize(aIndex+1);
        		array_[aIndex] = aValue;
        		return(size_);
    		}

    		if((size_ + 1) >= capacity_) {
        		int newCapacity;
        		bool success;
        		success = computeNewCapacity(size_ + 1,newCapacity);
        		if(!success) return(size_);
        		success = ensureCapacity(newCapacity);
        		if(!success) return(size_);
    		}

    		int i;
    		for(i = size_; i > aIndex; i--) {
        		array_[i] = array_[i-1];
    		}

    		array_[aIndex] = aValue;
    		size_++;

    		return(size_);
		}
		/**
 		 * @brief 	Remove a value from the array at a specified index.
 		 *
 		 * @details This method is relatively computationally costly since many of the array
 		 * 		elements may need to be shifted.
 		 *
 		 * @param aIndex Index of the value to remove.  All elements from aIndex to
 		 * 		the end of the array are shifted one place toward the beginning of
 		 * 		the array.  If aIndex is less than 0 or greater than or equal to the
 		 * 		current size of the array, no element is removed.
 		 * @return Size of the array after the removal.
 		*/
		int remove(int aIndex)
		{
    		if(aIndex < 0) {
        		std::cout << "Array.remove: ERR- aIndex was less than 0.\n";
        		return(size_);
    		}
    		if(aIndex >= size_) {
        		std::cout << "Array.remove: ERR- aIndex was greater than or equal the ";
        		std::cout << "size of the array.\n";
        		return(size_);
    		}

    		int i;
    		size_--;
    		for(i = aIndex; i < size_; i++) {
        		array_[i] = array_[i+1];
    		}
    		array_[size_] = defaultValue_;

    		return(size_);
		}
		/**
 		 * @brief Set the value at a specified index.
 		 *
 		 * @param aIndex Index of the array element to be set.  It is permissible
 		 * 		for aIndex to be past the current end of the array- the capacity will
 		 * 		be increased if necessary.  Values between the current end of the array
 		 * 		and aIndex are not initialized.
 		 * @param aValue Value.
 		*/
		void set(int aIndex,const T &aValue)
		{
    		if(aIndex < 0) return;

    		bool success = false;
    		if((aIndex+2) >= capacity_) {
        		int newCapacity;
        		success = computeNewCapacity(aIndex+2, newCapacity);
        		if(!success) return;
        		success = ensureCapacity(newCapacity);
        		if(!success) return;
    		}
    		array_[aIndex] = aValue;
    		if(aIndex >= size_)  size_ = aIndex + 1;
		}
		/**
 		 * @brief Get a pointer to the low-level array.
 		 *
 		 * @return Vecder to the low-level array.
 		*/
		T* get()
		{
    		return(array_);
		}
		/**
 		 * @brief Get a pointer to the low-level array.
 		 *
 		 * @return Vecder to the low-level array.
 		*/

		const T* get() const
		{
    		return(array_);
		}
		/**
 		 * @brief Get a const reference to the value at a specified array index.
 		 *
 		 * @details If the index is negative or passed the end of the array, an exception
 		 * is thrown.
 		 *
 		 * For faster execution, the array elements can be accessed through the
 		 * overloaded operator[], which does no bounds checking.
 		 *
 		 * @param aIndex Index of the desired array element.
 		 * @return const reference to the array element.
 		 * @throws Exception if (aIndex<0)||(aIndex>=size_).
 		 * @see operator[].
 		 */
		const T& get(int aIndex) const
		{
    		if((aIndex < 0) || (aIndex >= size_)) {
        		std::stringstream msg;
            	msg << "Array index out of bounds. " << ".";
            	throw (msg.str(),__FILE__,__LINE__);
    		}
    		return(array_[aIndex]);
		}
		/**
 		 * Get the last value in the array.
 		 *
 		 * @return Last value in the array.
 		 * @throws Exception if the array is empty.
 		*/
		const T& getLast() const
		{
    		if(size_ <= 0) {
        		std::stringstream msg;
            	msg << "Array is empty. " << ".";
            	throw (msg.str(),__FILE__,__LINE__);
    		}
    		return(array_[size_ - 1]);
		}
		/**
 		 * Get writable reference to last value in the array.
 		 *
 		 * @return writable reference to Last value in the array.
 		 * @throws Exception if the array is empty.
 		*/
		T& updLast() const
		{
    		if(size_ <= 0) {
        		std::stringstream msg;
            	msg << "Array is empty. " << ".";
            	throw (msg.str(),__FILE__,__LINE__);
    		}
    		return(array_[size_ - 1]);
		}
		/**
 		 * Linear search for an element matching a given value.
 		 *
 		 * @param aValue Value to which the array elements are compared.
 		 * @return Index of the array element matching aValue. If there is more than
 		 * one such elements with the same value the index of the first of these elements
 		 * is returned.  If no match is found, -1 is returned.
 		*/
		int findIndex(const T &aValue) const
		{
    		for(int i = 0; i < size_; i++) if(array_[i] == aValue) return i;
    		return -1;
		}

		/**
 		 * Linear search in reverse for an element matching a given value.
 		 *
 		 * @param aValue Value to which the array elements are compared.
 		 * @return Index of the array element matching aValue. If there is more than
 		 * one such elements with the same value the index of the last of these elements
 		 * is returned.  If no match is found, -1 is returned.
 		*/
		int rfindIndex(const T &aValue) const
		{
    		for(int i=size_ - 1; i >= 0; i--) if(array_[i]==aValue) return i;
    		return -1;
		}
	}; /** ended of class Array. */

}
#endif //ARRAY_H