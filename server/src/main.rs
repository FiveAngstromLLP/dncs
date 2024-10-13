use leptos::*;

#[component]
fn App() -> impl IntoView {
    let (count, set_count) = create_signal(0);

    view! {
        <div>
            <h1>"Welcome to Leptos!"</h1>
            <button
                on:click=move |_| set_count.update(|count| *count += 1)
            >
                "Click me: " {count}
            </button>
        </div>
    }
}

fn main() {
    mount_to_body(App)
}
