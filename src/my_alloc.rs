use std::alloc::{alloc, Layout, LayoutError};

struct MyPtr<T>{
    ptr: *mut T,
    elem_size: usize,
}

impl<T: Clone> MyPtr<T> {
    fn fill(&mut self, elem: T) {
        for i in 0..self.elem_size {
            unsafe {
                *self.ptr.wrapping_add(i) = elem.clone();
            }
        }
    }
}

fn my_alloc<T: Sized>(elem_size: usize) -> Result<MyPtr<T>, LayoutError> {
    let layout = Layout::from_size_align(elem_size * size_of::<T>(), size_of::<T>())?;
    let ptr = unsafe { alloc(layout) } as *mut T;
    Ok(MyPtr{ ptr, elem_size })
}