let keys = vec!["1fb13e", "1fb13d", "287146", "1fb13f", "1fb13a"];
let =  Vec<(i16, i16, i16)> = vec![(9435, 12302, 9174), (9408, 12297, 9166), (9959, 17480, 9016), (9463, 12307, 9186), (9331, 12296, 9164)];
let combined = keys.drain(..).zip(values.drain(..)).collect::<HashMap<_, (i16, i16, i16)>>();

