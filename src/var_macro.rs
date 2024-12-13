macro_rules! var {
    (mut $first_name: ident = $first_value: expr $(,$name: ident = $value:expr)*) => {
        let mut $first_name = $first_value;
        $(
            let mut $name = $value;
        )*
    };

    ($first_name: ident = $first_value: expr $(,$name: ident = $value:expr)*) => {
        let $first_name = $first_value;
        $(
            let $name = $value;
        )*
    };

    (mut $first_name: ident $(,$name:ident)* = $value:expr) => {
        let mut $first_name = $value;
        $(
            let mut $name = $value;
        )*
    };

    ($first_name: ident $(,$name:ident)* = $value:expr) => {
        let $first_name = $value;
        $(
            let $name = $value;
        )*
    };
}

macro_rules! te {
    ($ex:expr => $yes:expr ; $no:expr) => {
        if $ex { $yes } else { $no }
    }
}