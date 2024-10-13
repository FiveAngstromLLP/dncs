use leptos::*;

#[component]
fn App() -> impl IntoView {
    view! {
        <h1>"Hello, World!"</h1>
    }
}

fn main() {
    mount_to_body(App);
}
